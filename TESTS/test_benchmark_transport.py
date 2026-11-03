#!/usr/bin/env python3
"""Check benchmark diagnostic modes and rejection gates with synthetic outputs."""
import argparse
import copy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import f90nml

import benchmark_transport as benchmark


class BenchmarkDiagnostics(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.case = Path(self.temp.name)
        self.args = argparse.Namespace(collisions="lorentz", time_step=.25)
        self.outputs = {
            "transport_metadata.dat": ".1 2700 .9 1.1 3e9 1e5 .25 4 2 2 5\n",
            "n_lost_particles.dat": "1\n",
            "loss_summary.dat": ".1 1 2\n",
            "lost_particle_events.dat": ".1 1 1 2 .5 1 0 0\n",
            "nustar_diffcoef_std.dat": ".1 1 .1\n",
            "mean_square_displacement.dat": ".25 .1\n.5 .2\n.75 .3\n1 .4\n0 0\n",
            "mean_square_uncertainty.dat": ".25 .01\n.5 .02\n.75 .03\n1 .04\n0 0\n",
            "particle_histories.dat": "".join(
                f".1 {marker} {step} {step * .25} .5\n"
                for marker, steps in ((1, range(2)), (2, range(5))) for step in steps)}
        self.record = {
            "errors": [], "particles": 2, "backend": "cpu-batch", "returncode": 0,
            "counters": {"attempted": 2, "completed": 1, "physically_lost": 1,
                         "numerically_failed": 0, "resolved_steps": 4,
                         "orbit_advances": 6, "cell_crossings": 8, "replay_cells": 0},
            "profile": {name: 0. for name in ("setup", "upload", "kernel", "download", "replay", "wall")}}

    def write_outputs(self):
        for name, contents in self.outputs.items():
            (self.case / name).write_text(contents)

    def errors(self):
        record = copy.deepcopy(self.record)
        benchmark.validate(self.case, self.args, record)
        return record["errors"]

    def test_flag_and_prepared_inputs(self):
        argv = ["benchmark_transport.py", "--binary", "binary", "--fixture", str(self.case),
                "--output", "output", "--collisions", "lorentz"]
        for disabled in (False, True):
            with self.subTest(disabled=disabled), patch.object(
                    sys, "argv", argv + (["--no-particle-histories"] if disabled else [])):
                args = benchmark.arguments()
                self.assertEqual(args.no_particle_histories, disabled)
                for name in benchmark.INPUTS:
                    (self.case / name).write_text("")
                for name, group in (("gorilla.inp", "gorillanml"),
                                    ("tetra_grid.inp", "tetra_grid_nml"),
                                    ("gorilla_applets.inp", "gorilla_applets_nml"),
                                    ("mono_energetic_transp_coef.inp", "transpcoefnml")):
                    (self.case / name).write_text(f"&{group}\n/\n")
                (self.case / "MHD_EQUILIBRIA").mkdir(exist_ok=True)
                for name in benchmark.EQUILIBRIA:
                    (self.case / "MHD_EQUILIBRIA" / name).write_text("")
                prepared = self.case / str(disabled)
                benchmark.prepare(prepared, args, 2, "cpu-batch")
                mono = f90nml.read(prepared / "mono_energetic_transp_coef.inp")["transpcoefnml"]
                self.assertEqual(mono["boole_write_particle_histories"], not disabled)
                self.assertEqual(mono["idiffcoef_output"], 3)
                self.assertTrue(mono["boole_psi_mat"])

    def test_default_requires_complete_histories(self):
        self.write_outputs()
        self.assertFalse(self.errors())
        path = self.case / "particle_histories.dat"
        path.write_text(self.outputs[path.name].split("\n", 1)[1])
        self.assertTrue(self.errors())
        path.unlink()
        self.assertTrue(self.errors())
        self.args.no_particle_histories = True
        self.assertFalse(self.errors())
        self.assertNotIn(path.name, benchmark.output_files(self.args))

    def test_corrupt_outputs_rejected_in_both_modes(self):
        corruptions = {
            "transport_metadata.dat": ["1 2\n", ".1 2700 .9 1.1 3e9 1e5 .25 4 2 1 5\n",
                                       ".1 2700 .9 1.1 3e9 1e5 .5 4 2 2 5\n"],
            "lost_particle_events.dat": [".1 1 1 2\n", ".1 3 1 2 .5 1 0 0\n",
                                         ".1 1 1.5 2 .5 1 0 0\n", ".1 1 1 5 1.25 1 0 0\n"],
            "nustar_diffcoef_std.dat": [".1 NaN .1\n"],
            "mean_square_displacement.dat": [".25 .1\n0 0\n"],
            "mean_square_uncertainty.dat": [".25 NaN\n"]}
        for disabled in (False, True):
            self.args.no_particle_histories = disabled
            for name, variants in corruptions.items():
                for contents in variants:
                    with self.subTest(disabled=disabled, name=name, contents=contents):
                        self.write_outputs()
                        (self.case / name).write_text(contents)
                        self.assertTrue(self.errors())

    def test_loss_identifiers_compared_exactly_in_both_modes(self):
        reference = self.case / "reference"
        reference.mkdir()
        self.write_outputs()
        for name, contents in self.outputs.items():
            (reference / name).write_text(contents)
        for disabled in (False, True):
            self.args.no_particle_histories = disabled
            for column in (1, 2, 3):
                with self.subTest(disabled=disabled, column=column):
                    row = self.outputs["lost_particle_events.dat"].split()
                    row[column] = str(float(row[column]) + 1.e-10)
                    (self.case / "lost_particle_events.dat").write_text(" ".join(row) + "\n")
                    errors, _ = benchmark.compare(self.case, reference, self.args)
                    self.assertTrue(errors)

    def test_first_sample_and_repeat_statistics(self):
        for mode in ("full-particle-histories", "moments-and-losses"):
            records = [{"particles": 2, "backend": "cpu-batch", "threads": 1,
                        "diagnostic_mode": mode, "phase": phase, "wall_seconds": wall,
                        "valid": True, "speedup_eligible": True}
                       for phase, wall in (("cold_process", 100.), ("repeat_process", 2.),
                                           ("repeat_process", 4.))]
            summary = benchmark.summarize(records)[0]
            self.assertEqual(summary["diagnostic_mode"], mode)
            self.assertEqual(summary["repeats"], 2)
            self.assertEqual(summary["median_wall_seconds"], 3.)
            self.assertEqual(summary["min_wall_seconds"], 2.)
            self.assertEqual(summary["max_wall_seconds"], 4.)
            records[0]["valid"] = False
            summary = benchmark.summarize(records)[0]
            self.assertIsNone(summary["median_wall_seconds"])
            self.assertIsNone(summary["speedup_vs_best_measured_cpu"])


if __name__ == "__main__":
    unittest.main()
