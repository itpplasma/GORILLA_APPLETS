# TODO: GORILLA Kramers-Moyal D11 Benchmark

Branch: `transport-km-benchmark`
Benchmark target: GORILLA vs NEO-2-QL vs NEOART on AUG 39461 neon discharge

## Done

- [x] KM unit tests: 1D random walk recovery, fit_transport_coefficients validation
- [x] KM module: `kramers_moyal_transport_mod.f90` wrapping existing displacement infrastructure
- [x] Updated CLAUDE.md with KM approach, archived LM references
- [x] Branch created and pushed, stacked on `transport-tracing-guard` (PR #22)

## Remaining Work

### 1. Wire KM module into i_option or standalone driver

Add a new i_option (e.g. 16) or a standalone benchmark subroutine that:
- Reads equilibrium from EQDSK/Boozer (AUG 39461: `$DATA/AUG/BOOZER/39461/39461.bc`)
- Initializes GORILLA with the tokamak grid
- Calls `calc_km_d11_profile` at the same 20 flux surfaces as NEO-2-QL (`boozer_s` from HDF5)
- Writes D11(s) to CSV

Files:
- `SRC/gorilla_applets_main.f90`: add case for new i_option
- `SRC/kramers_moyal_transport_mod.f90`: may need to accept surface s-values directly (not just indices)

### 2. Geometry conversion: D11_s to D11_AX

The KM module outputs D11 in s-coordinate units. NEO-2 reports D11_AX in cm^2/s.

Conversion: `D11_AX = D11_s / <nabla_s>^2`

where `<nabla_s>` = `av_nabla_stor` from NEO-2 output or computed from GORILLA geometry.

Formula derivation in `$DATA/AUG/NEO-2/neo2_neoart_benchmark_tools/doc/README.md`.

Need to either:
- Read `av_nabla_stor` from the NEO-2 HDF5 output, or
- Compute `<nabla_s>` from the GORILLA tetrahedral mesh geometry

Files:
- `SRC/kramers_moyal_transport_mod.f90`: add conversion routine
- Test: verify conversion with known S, V, R_axis values

### 3. Benchmark run on AUG 39461

Run GORILLA KM D11 on the 20 NEO-2 flux surfaces for AUG 39461.

Input data:
- Equilibrium: `$DATA/AUG/EQDSK/39461/eqdsk_39461_5.38s`
- Boozer: `$DATA/AUG/BOOZER/39461/39461.bc`
- Grid: need to determine appropriate n1/n2/n3 for this equilibrium

Parameters:
- Species: electrons (tracer_species=1)
- Collision operator: Lorentz (collision_operator=4, mono-energetic)
- Energy: match NEO-2 (check `T_spec` in HDF5, likely 3.5 keV for electrons)
- v_E: 0 (no ExB drift initially)
- Particles: 5000-40000 per surface (start with 5000, increase for convergence)
- Tracing time: 2x collision time (auto-estimated from nu_star)

NEO-2-QL reference:
- `$DATA/AUG/NEO-2/39461/neoart_benchmark_neon_discharge/neon_discharge_out.h5`
- 20 surfaces: boozer_s = 0.1 to 0.98
- D11_AX shape (20, 25) -- 20 surfaces, 25 species pairs
- av_nabla_stor available for conversion

NEOART reference:
- `$DATA/AUG/NEO-ART/39461/neoart_benchmark_neon_discharge/neon_discharge_neoart.nc`

Output directory: `$DATA/AUG/GORILLA/39461/neoart_benchmark_neon_discharge/`

Script: `scripts/run_km_benchmark_39461.sh`

### 4. Comparison plots

Generate plots comparing GORILLA D11 against NEO-2-QL and NEOART:

1. **D11(rho_pol) radial profile**: all three codes on same axes
2. **Convergence study**: D11 vs N_particles at one surface (e.g. s=0.5)
3. **A(s) drift coefficient**: verify A ~= 0 for Lorentz with v_E=0

Use existing comparison tools at `$DATA/AUG/NEO-2/neo2_neoart_benchmark_tools/` or write new plotting script.

Upload to Litterbox, link in PR.

### 5. Cleanup stale data

Remove (user must do manually):
```
rm -rf ~/data/AUG/GORILLA/30835
rm -rf ~/data/AUG/NEO-2/30835
rm -rf ~/data/AUG/COMPARISONS/30835
```

### 6. Update $DATA/AUG/README.md

Add GORILLA benchmark project to the project list:
```
- Benchmark of mono-energetic D11 between GORILLA (Kramers-Moyal) and NEO-2-QL
  for AUG shot 39461

    AUG/
    ├── GORILLA/
    │   └── 39461/neoart_benchmark_neon_discharge/
    ├── NEO-2/
    │   └── 39461/neoart_benchmark_neon_discharge/
    └── NEO-ART/
        └── 39461/neoart_benchmark_neon_discharge/
```

### 7. Documentation in $DATA benchmark directories

- `$DATA/AUG/GORILLA/39461/neoart_benchmark_neon_discharge/README.md`: document GORILLA setup, parameters, how to reproduce
- Reference the NEO-2/NEOART benchmark tools documentation

### 8. PR

Create PR for `transport-km-benchmark` with:
- Summary of KM approach vs abandoned LM approach
- D11 comparison plots (Litterbox links)
- Test results
- Link to benchmark data location

## Physics Notes

- D11 = <(ds)^2>/(2dt) at each flux surface (Kramers-Moyal 2nd coefficient)
- A = <ds>/dt at each flux surface (1st coefficient, should be ~0 for Lorentz/v_E=0)
- Conversion to physical units: D11_AX = D11_s / <nabla_s>^2
- NEO-2 uses NEO-2-QL (field-line integration, tokamak version) for this benchmark
- NEOART is a fluid code (banana-plateau + Pfirsch-Schluter), expected to match NEO-2 in high-collisionality edge
- GORILLA uses particle orbit tracing with Lorentz collision operator (mono-energetic)

## Key Files

| File | Purpose |
|------|---------|
| `SRC/kramers_moyal_transport_mod.f90` | KM D11 profile estimation |
| `SRC/transport_benchmark_utils_mod.f90` | fit_transport_coefficients (KM fitting) |
| `SRC/utils_self_consistent_ef_mod.f90:1254-1422` | Original D11 estimation loop |
| `SRC/utils_self_consistent_ef_mod.f90:377-399` | Displacement accumulation in pusher |
| `SRC/gorilla_applets_types_mod.f90:207-234` | delta_s_delta_s_squared_t type |
| `TESTS/test_kramers_moyal.f90` | KM unit tests |
| `$DATA/AUG/NEO-2/neo2_neoart_benchmark_tools/doc/README.md` | Conversion formulas |
| `$DATA/AUG/NEO-2/39461/.../neon_discharge_out.h5` | NEO-2-QL reference D11 |
| `$DATA/AUG/NEO-ART/39461/.../neon_discharge_neoart.nc` | NEOART reference D11 |
