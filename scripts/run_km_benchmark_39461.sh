#!/usr/bin/env bash
set -euo pipefail

# GORILLA KM D11 benchmark for AUG shot 39461 (neon discharge)
# Compares mono-energetic electron D11 against NEO-2-QL

DATA="${DATA:-${HOME}/data}"
repo_root="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
gorilla_root="${repo_root%/GORILLA_APPLETS}/GORILLA"
case_dir="${1:-${DATA}/AUG/GORILLA/39461/neoart_benchmark_neon_discharge}"
n_particles="${N_PARTICLES:-5000}"
grid_n1="${GRID_N1:-24}"
grid_n2="${GRID_N2:-20}"
grid_n3="${GRID_N3:-24}"

eqdsk="${DATA}/AUG/EQDSK/39461/eqdsk_39461_5.38s"
convex_wall="${DATA}/AUG/BOOZER/39461/convexwall.dat"
neo2_h5="${DATA}/AUG/NEO-2/39461/neoart_benchmark_neon_discharge/neon_discharge_out.h5"
neoart_nc="${DATA}/AUG/NEO-ART/39461/neoart_benchmark_neon_discharge/neon_discharge_neoart.nc"

for f in "$eqdsk" "$convex_wall" "$neo2_h5"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: required file not found: $f"
        exit 1
    fi
done

if [[ ! -f "${gorilla_root}/954/F90/Src/Polynomial234RootSolvers.f90" ]]; then
    "${gorilla_root}/fetch_polynomial234roots.sh"
fi

make -C "${gorilla_root}"
make -C "${repo_root}"

# Step 1: Precompute fluxtube volumes
fluxtv_dir="${case_dir}/fluxtv_precomp"
mkdir -p "${fluxtv_dir}"

ln -snf "$(realpath --relative-to "${fluxtv_dir}" "${repo_root}/BUILD/gorilla_applets_main.x")" \
    "${fluxtv_dir}/gorilla_applets_main.x"
ln -snf "$(realpath --relative-to "${fluxtv_dir}" "${gorilla_root}/INPUT/preload_for_SYNCH.inp")" \
    "${fluxtv_dir}/preload_for_SYNCH.inp"
ln -snf "$(realpath --relative-to "${fluxtv_dir}" "${repo_root}/INPUT/seed.inp")" \
    "${fluxtv_dir}/seed.inp"

cat > "${fluxtv_dir}/field_divB0.inp" <<'EOF'
0                                 ipert
1                                 iequil
1.00                              ampl
72                                ntor
0.99                              cutoff
4                                 icftype
'see_tetra_grid_inp'              gfile
'dummy'                           pfile
'see_tetra_grid_inp'              convexfile
'DATA/ASDEX/FLUXDATA'             fluxdatapath
0                                 window size for filtering of psi array over R
0                                 window size for filtering of psi array over Z
0                                 iaxieq
EOF

cat > "${fluxtv_dir}/tetra_grid.inp" <<EOF
&TETRA_GRID_NML
  n1 = ${grid_n1} ,
  n2 = ${grid_n2} ,
  n3 = ${grid_n3} ,
  grid_kind = 2 ,
  g_file_filename = '${eqdsk}' ,
  convex_wall_filename = '${convex_wall}' ,
  netcdf_filename = 'wout_dummy.nc' ,
  knots_SOLEDGE3X_EIRENE_filename = 'dummy' ,
  triangles_SOLEDGE3X_EIRENE_filename = 'dummy' ,
  boole_n_field_periods = .true. ,
  n_field_periods_manual = 1 ,
  sfc_s_min = 1.d-4 ,
  theta_geom_flux = 1 ,
  theta0_at_xpoint = .true. ,
  boole_write_mesh_obj = .false. ,
  filename_mesh_rphiz = 'mesh_rphiz.obj' ,
  filename_mesh_sthetaphi = 'mesh_sthetaphi.obj'
/
EOF

cat > "${fluxtv_dir}/gorilla.inp" <<'EOF'
&GORILLANML
  eps_Phi = 0.0d0 ,
  coord_system = 2 ,
  ispecies = 1 ,
  boole_periodic_relocation = .true. ,
  ipusher = 2 ,
  boole_pusher_ode45 = .false. ,
  rel_err_ode45 = 1.E-10 ,
  boole_dt_dtau = .true. ,
  boole_newton_precalc = .false. ,
  poly_order = 2 ,
  i_precomp = 0 ,
  boole_guess = .true. ,
  i_time_tracing_option = 1 ,
  handover_processing_kind = 1 ,
  boole_time_hamiltonian = .true. ,
  boole_gyrophase = .true. ,
  boole_vpar_int = .true. ,
  boole_vpar2_int = .true. ,
  boole_adaptive_time_steps = .false. ,
  desired_delta_energy = 1.E-10 ,
  max_n_intermediate_steps = 10000 ,
  boole_grid_for_find_tetra = .true. ,
  a_factor = 1 ,
  b_factor = 0 ,
  c_factor = 0 ,
  boole_strong_electric_field = .false. ,
  boole_save_electric = .false. ,
  filename_electric_field = 'electric_field.dat' ,
  filename_electric_drift = 'electric_drift.dat' ,
  boole_pert_from_mephit = .false.
/
EOF

cat > "${fluxtv_dir}/gorilla_applets.inp" <<'EOF'
&GORILLA_APPLETS_NML
  i_option = 1 ,
  filename_fluxtv_precomp = 'fluxtubevolume.dat' ,
  start_pos_x1 = 0.5d0 ,
  start_pos_x2 = 0.0d0 ,
  start_pos_x3 = 0.0d0 ,
  t_step_fluxtv = 1.d-4 ,
  nt_steps_fluxtv = 100000 ,
  energy_eV_fluxtv = 5.0d2 ,
  filename_fluxtv_load = 'fluxtubevolume.dat'
/
EOF

printf 'Precomputing fluxtube volumes for AUG 39461 ...\n'
(
    cd "${fluxtv_dir}"
    OMP_NUM_THREADS="${OMP_THREADS:-$(nproc)}" ./gorilla_applets_main.x > run.log 2>&1
)
if [[ ! -f "${fluxtv_dir}/fluxtubevolume.dat" ]]; then
    echo "ERROR: fluxtube volume precomputation failed"
    tail -20 "${fluxtv_dir}/run.log"
    exit 1
fi
printf 'Fluxtube volume precomputation done.\n'

# Step 2: Run KM benchmark
run_dir="${case_dir}/km_benchmark"
mkdir -p "${run_dir}"

ln -snf "$(realpath --relative-to "${run_dir}" "${repo_root}/BUILD/gorilla_applets_main.x")" \
    "${run_dir}/gorilla_applets_main.x"
ln -snf "$(realpath --relative-to "${run_dir}" "${repo_root}/INPUT/seed.inp")" \
    "${run_dir}/seed.inp"
ln -snf "$(realpath --relative-to "${run_dir}" "${gorilla_root}/INPUT/preload_for_SYNCH.inp")" \
    "${run_dir}/preload_for_SYNCH.inp"

cp "${fluxtv_dir}/field_divB0.inp" "${run_dir}/"
cp "${fluxtv_dir}/tetra_grid.inp" "${run_dir}/"
cp "${fluxtv_dir}/gorilla.inp" "${run_dir}/"
cp "${fluxtv_dir}/fluxtubevolume.dat" "${run_dir}/"

cat > "${run_dir}/gorilla_applets.inp" <<'EOF'
&GORILLA_APPLETS_NML
  i_option = 14 ,
  filename_fluxtv_precomp = 'fluxtubevolume.dat' ,
  start_pos_x1 = 0.5d0 ,
  start_pos_x2 = 0.0d0 ,
  start_pos_x3 = 0.0d0 ,
  t_step_fluxtv = 1.d-4 ,
  nt_steps_fluxtv = 100000 ,
  energy_eV_fluxtv = 5.0d2 ,
  filename_fluxtv_load = 'fluxtubevolume.dat'
/
EOF

cat > "${run_dir}/km_benchmark.inp" <<EOF
&km_benchmark_nml
  n_particles = ${n_particles} ,
  tracer_species = 1 ,
  collision_operator = 4 ,
  boole_precalc_collisions = .false. ,
  temperature_eV = 5.0d2 ,
  total_time = 0.0d0 ,
  v_E = 0.0d0 ,
  i_integrator_type = 1 ,
  n_background_species = 5 ,
  background_mass_amu = 2.014d0, 20.18d0, 20.18d0, 20.18d0, 5.486d-4 ,
  background_charge = 1.0d0, 10.0d0, 9.0d0, 8.0d0, -1.0d0 ,
  background_density = 7.1d13, 3.0d11, 1.8d11, 9.5d10, 7.7d13 ,
  background_temperature = 5.0d2, 5.0d2, 5.0d2, 5.0d2, 5.0d2 ,
  n_surfaces = 20 ,
  surface_s_values = 0.10000000d0, 0.14631579d0, 0.19263158d0, 0.23894737d0, 0.28526316d0,
                     0.33157895d0, 0.37789474d0, 0.42421053d0, 0.47052632d0, 0.51684211d0,
                     0.56315789d0, 0.60947368d0, 0.65578947d0, 0.70210526d0, 0.74842105d0,
                     0.79473684d0, 0.84105263d0, 0.88736842d0, 0.93368421d0, 0.98000000d0 ,
  filename_output = 'km_d11_profile.csv' ,
  collision_profile_file = 'neo2_collision_profile_39461.dat'
/
EOF

# Generate collision profile from NEO-2 data
if [[ -f "${repo_root}/.venv/bin/python3" ]]; then
    py="${repo_root}/.venv/bin/python3"
else
    py=python3
fi
"${py}" -c "
import h5py, numpy as np
e=4.8032e-10; m=9.1094e-28; ev=1.6022e-12
f=h5py.File('${neo2_h5}','r')
bs=f['boozer_s'][:]; ns=f['n_spec'][:]; cp=f['collpar_spec'][:,0]; f.close()
with open('${run_dir}/neo2_collision_profile_39461.dat','w') as o:
    o.write('# NEO-2 collision profile AUG 39461\n')
    o.write('# Species: D+(1), Ne10+(10), Ne9+(9), Ne8+(8), e-(-1)\n')
    o.write(f'{len(bs)} 5\n')
    o.write('# s  n_D  n_Ne10  n_Ne9  n_Ne8  n_e  E_eV\n')
    for i in range(len(bs)):
        z2n=ns[i,1]+ns[i,2]*100+ns[i,3]*81+ns[i,4]*64
        v=(2*np.pi*e**4*z2n*14/(m**2*cp[i]))**0.25
        E=0.5*m*v**2/ev
        o.write(f'{bs[i]:.8f}  {ns[i,1]:.6e}  {ns[i,2]:.6e}  {ns[i,3]:.6e}  {ns[i,4]:.6e}  {ns[i,0]:.6e}  {E:.1f}\n')
"

printf 'Running KM D11 benchmark on AUG 39461 (%s particles) ...\n' "${n_particles}"
(
    cd "${run_dir}"
    OMP_NUM_THREADS="${OMP_THREADS:-$(nproc)}" ./gorilla_applets_main.x > run.log 2>&1
)

if [[ ! -f "${run_dir}/km_d11_profile.csv" ]]; then
    echo "ERROR: KM benchmark failed"
    tail -30 "${run_dir}/run.log"
    exit 1
fi

printf 'KM benchmark done. Results in %s\n' "${run_dir}/km_d11_profile.csv"
cat "${run_dir}/km_d11_profile.csv"

# Step 3: Generate comparison plot
plot_args=(
    "${run_dir}/km_d11_profile.csv"
    "${neo2_h5}"
    --output "${case_dir}/gorilla_d11_comparison.png"
)
if [[ -f "${neoart_nc}" ]]; then
    plot_args+=(--neoart-nc "${neoart_nc}")
fi

printf 'Generating comparison plot ...\n'
if [[ -f "${repo_root}/.venv/bin/python3" ]]; then
    "${repo_root}/.venv/bin/python3" "${repo_root}/scripts/convert_and_compare_d11.py" "${plot_args[@]}"
else
    python3 "${repo_root}/scripts/convert_and_compare_d11.py" "${plot_args[@]}"
fi
printf 'Done. Plot: %s\n' "${case_dir}/gorilla_d11_comparison.png"
