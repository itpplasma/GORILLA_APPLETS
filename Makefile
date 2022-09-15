FC = gfortran
#OPTS= -M OBJS --chk a,e,s,u,x --trace --trap -g
#OPTS= -M OBJS -O 
#OPTS= -module OBJS
#OPTS= -J OBJS -O
#vOPTS= -J OBJS -O0 -g -fbounds-check
#OPTS= -J OBJS  -g -fbacktrace -ffpe-trap=zero,overflow,invalid  -fbounds-check
#OPTS = -J OBJS -Wall -pedantic
#OPTS = -J OBJS -Wuninitialized
#OPTS = -J OBJS -g

# The developer repo option
OPTS= -g -J OBJS -O3 -fopenmp 
# The option below shows in which line the error occurs
#OPTS= -J OBJS  -g -fbacktrace -ffpe-trap=zero,overflow,invalid  -fbounds-check -fcheck=all -fopenmp
#OPTS= -J OBJS  -ffpe-trap=invalid,zero,overflow -g -fopenmp
# Official Core module options
#OPTS ?= -J OBJS  -g -fbacktrace -ffpe-trap=zero,overflow,invalid  -fbounds-check -fopenmp

UNAME_S	:= $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	NCINC ?= -I/opt/local/include
	NCLIB ?= -L/opt/local/lib -lnetcdff -lnetcdf -llapack
else
	NCINC ?= -I/usr/include
	NCLIB ?= -lnetcdff -lnetcdf -llapack
endif

SOURCES = SetWorkingPrecision.f90\
	Polynomial234RootSolvers.f90 \
	llsq.f90 \
	constants_mod.f90 \
	tetra_grid_settings_mod.f90 \
	gorilla_settings_mod.f90 \
	gorilla_applets_settings_mod.f90 \
	mono_energetic_transp_coef_settings_mod.f90 \
	various_functions_mod.f90 \
	gorilla_diag_mod.f90 \
	canonical_coordinates_mod.f90 \
	nctools_module.f90 \
	rkf45.f90 \
	odeint_rkf45.f90 \
	runge_kutta_mod.f90 \
	magfie.f90 \
	chamb_m.f90 \
	vmecinm_m.f90 \
	spl_three_to_five_mod.f90 \
	spline_vmec_data.f90 \
	new_vmec_allocation_stuff.f90 \
	binsrc.f90 \
	field_divB0.f90 \
	scaling_r_theta.f90\
	field_line_integration_for_SYNCH.f90 \
	preload_for_SYNCH.f90 \
	plag_coeff.f90 \
	magdata_in_symfluxcoord.f90 \
	points_2d.f90\
	circular_mesh.f90\
	circular_mesh_SOLEDGE3X_EIRENE.f90\
	tetra_grid_mod.f90 \
	make_grid_rect.f90 \
	bdivfree.f90 \
	tetra_physics_mod.f90 \
	tetra_physics_poly_precomp_mod.f90 \
	sub_alpha_lifetime_can.f90 \
	differentiate.f90 \
	spline5_RZ.f90 \
	supporting_functions_mod.f90 \
	pusher_tetra_func_mod.f90 \
	pusher_tetra_poly.f90 \
	pusher_tetra_rk.f90 \
	get_canonical_coordinates.f90 \
	orbit_timestep_gorilla.f90 \
	gorilla_plot_mod.f90 \
	fluxtv_pusher_mod.f90 \
	fluxtv_mod.f90 \
	flux_deviation_mod.f90 \
	gorilla_applets_sub_mod.f90 \
	alpha_lifetime_gorilla_mod.f90 \
	direct_vmec_integrator_mod.f90 \
	poincare_invariances_mod.f90 \
	boltzmann_mod.f90 \
	total_dwell_times_mod.f90 \
	gorilla_applets_main.f90

OBJS = $(patsubst %.f90,OBJS/%.o,$(SOURCES))


gorilla_applets_main.x: $(OBJS_CONTRIB) $(OBJS)
	$(FC) $(OPTS) -o $@ $^ $(NCLIB)

OBJS/%.o: SRC/contrib/%.f90
	$(FC) $(OPTS) -c $^ -o $@ $(NCINC)

OBJS/%.o: SRC/%.f90
	$(FC) $(OPTS) -c $^ -o $@ $(NCINC)

OBJS/%.o: SRC_CORE/contrib/%.f90
	$(FC) $(OPTS) -c $^ -o $@ $(NCINC)

OBJS/%.o: SRC_CORE/%.f90
	$(FC) $(OPTS) -c $^ -o $@ $(NCINC)


.PHONY: clean
clean:
	rm -f OBJS/*
	rm -f SRC/*.mod
	rm -f SRC_CORE/*.mod
	rm -f gorilla_applets_main.x
