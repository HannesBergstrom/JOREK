# --- General settings
jorekmodel="600"
description="Kinetic RE drift-surface equilibrium (re_kinetic_equilibrium) with model$jorekmodel: multi-class RE current on constant-P_phi surfaces, q-profile matching, and marker initialization from the converged equilibrium via kinetic_main (rep coupling scheme, full-orbit markers)."
options="with_vpar=.false. with_TiTe=.false. with_neutrals=.false. with_impurities=.false. with_refluid=.false."
mpitasks=2
num_threads=8
binaries="kinetic_main"
binaries_initial="jorek_model${jorekmodel}_1"
extra_restart="part_restart.h5 re_equilibrium.dat"
requiredfiles="input re_distribution.dat qprofile_target.dat"
extra_remote_files="part_restart.h5 re_equilibrium.dat"
particle_example="kinetic_main"
particle_example_dir="particles/examples"


# --- Compile the code for the test case
function compile_jorek () {
  if [ "$initialrun" == "yes" ]; then
    ./util/config.sh model=$jorekmodel n_tor=1 n_coord_tor=1 l_pol_domm=0 n_period=1 n_coord_period=1 n_plane=1 n_order=3 $options || exit 1
    make $compilopt $debugoptions jorek_model${jorekmodel}                                                                         || exit 1
    mv jorek_model${jorekmodel} jorek_model${jorekmodel}_1                                                                         || exit 1
    make cleanall                                                                                                                  || exit 1
  fi
  ./util/config.sh model=$jorekmodel n_tor=3 n_coord_tor=1 l_pol_domm=0 n_plane=4 n_period=1 n_coord_period=1 n_order=3 $options   || exit 1
  make $compilopt $debugoptions kinetic_main                                                                                       || exit 1
}


# --- Initial run only required when preparing or updating the test case
function initial_run () {
  # 1) drift-surface equilibrium phase: writes jorek00000.h5 and re_equilibrium.dat
  #    (re_eq_convergence.log documents the q-matching; must end converged)
  ${codedir}/util/setinput.sh input restart_particles=.f. restart=.f. nstep_n=0 tstep_n=1. linear_run=.t.    || exit 1
  export OMP_NUM_THREADS=$num_threads                                                                        || exit 1
  echo "setting OMP_NUM_THREADS=$num_threads, due to the requirements of the test"                           || exit 1
  $MPIRUN $mpitasks ./jorek_model${jorekmodel}_1 < input | tee logfile_initial                               || exit 1
  # 2) marker initialization from the equilibrium + a short coupled run:
  #    writes part_restart_init.h5 / part_restart.h5
  ${codedir}/util/setinput.sh input restart_particles=.f. restart=.t. nstep_n=10 tstep_n=0.01 linear_run=.t. || exit 1
  $MPIRUN $mpitasks ./kinetic_main < input | tee logfile_initial2
}


# --- Carry out the test case
function restart_run () {
  ${codedir}/util/setinput.sh input restart_particles=.t. restart=.t. nstep_n=2 tstep_n=0.5 nout=1 || exit 1
  export OMP_NUM_THREADS=$num_threads                                                              || exit 1
  echo "setting OMP_NUM_THREADS=$num_threads, due to the requirements of the test"                 || exit 1
  $MPIRUN $mpitasks ./kinetic_main < input | tee logfile                                           || exit 1
}


# --- Compare the results of the test case to the reference solution
function compare_results () {
  compare_results_generic 1.e-10 || exit 1
}
