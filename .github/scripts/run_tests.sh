#!/bin/bash
#
# runs a test example case
#

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

WORKDIR=`pwd`
dir=${TESTDIR}

# info
echo "work directory: $WORKDIR"
echo `date`
echo
echo "**********************************************************"
echo
echo "test directory: $dir"
echo
echo "**********************************************************"
echo

# bash function for checking seismogram output with reference solutions
my_test(){
  echo "testing seismograms:"
  ln -s $WORKDIR/utils/scripts/compare_seismogram_correlations.py
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/
  if [[ $? -ne 0 ]]; then exit 1; fi
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/ | grep min/max | cut -d \| -f 3 | awk '{print "correlation:",$1; if ($1 < 0.999 ){print $1,"failed"; exit 1;}else{ print $1,"good"; exit 0;}}'
  if [[ $? -ne 0 ]]; then exit 1; fi
}

my_kernel_test(){
  # kernel value test - checks rho/kappa/mu kernel value outputs
  echo "testing kernel values:"
  if [ ! -e REF_KERNEL/output_solver.txt ]; then echo "Please check if file REF_KERNEL/output_solver.txt exists..."; ls -alR ./; exit 1; fi
  if [ ! -e OUTPUT_FILES/output_solver.txt ]; then echo "Please check if file OUTPUT_FILES/output_solver.txt exists..."; ls -alR ./; exit 1; fi
  # gets reference expected kernel values from REF_KERNEL/ folder
  RHO=`grep -E 'maximum value of rho[[:space:]]+kernel' REF_KERNEL/output_solver.txt | cut -d = -f 2 | tr -d ' '`
  KAPPA=`grep -E 'maximum value of kappa[[:space:]]+kernel' REF_KERNEL/output_solver.txt | cut -d = -f 2 | tr -d ' '`
  MU=`grep -E 'maximum value of mu[[:space:]]+kernel' REF_KERNEL/output_solver.txt | cut -d = -f 2 | tr -d ' '`
  # need at least rho & kappa (for acoustic kernels)
  if [ "$RHO" == "" ] || [ "$KAPPA" == "" ]; then
    echo "  missing reference kernel values: RHO=$RHO KAPPA=$KAPPA MU=$MU"
    echo
    exit 1
  else
    echo "  reference kernel values: RHO=$RHO KAPPA=$KAPPA MU=$MU"
    echo
  fi
  # compares with test output
  # final test result
  PASSED=0
  # checks rho kernel value
  if [ "$RHO" != "" ]; then
    VAL=`grep -E 'maximum value of rho[[:space:]]+kernel' OUTPUT_FILES/output_solver.txt | cut -d = -f 2 | tr -d ' '`
    echo "kernel rho   : $VAL"
    echo "" | awk '{diff=ex-val;diff_abs=(diff >= 0)? diff:-diff;diff_rel=diff_abs/ex;print "  value: expected = "ex" gotten = "val" - difference absolute = "diff_abs" relative = "diff_rel; if (diff_rel>0.0001){print "  failed"; exit 1;}else{print "  good"; exit 0;} }' ex=$RHO val=$VAL
    if [[ $? -ne 0 ]]; then PASSED=1; fi
  fi
  # checks kappa kernel value
  if [ "$KAPPA" != "" ]; then
    VAL=`grep -E 'maximum value of kappa[[:space:]]+kernel' OUTPUT_FILES/output_solver.txt | cut -d = -f 2 | tr -d ' '`
    echo "kernel kappa : $VAL"
    echo "" | awk '{diff=ex-val;diff_abs=(diff >= 0)? diff:-diff;diff_rel=diff_abs/ex;print "  value: expected = "ex" gotten = "val" - difference absolute = "diff_abs" relative = "diff_rel; if (diff_rel>0.0001){print "  failed"; exit 1;}else{print "  good"; exit 0;} }' ex=$KAPPA val=$VAL
    if [[ $? -ne 0 ]]; then PASSED=1; fi
  fi
  # checks mu kernel value (if available for elastic kernel)
  if [ "$MU" != "" ]; then
    VAL=`grep -E 'maximum value of mu[[:space:]]+kernel' OUTPUT_FILES/output_solver.txt | cut -d = -f 2 | tr -d ' '`
    echo "kernel mu    : $VAL"
    echo "" | awk '{diff=ex-val;diff_abs=(diff >= 0)? diff:-diff;diff_rel=diff_abs/ex;print "  value: expected = "ex" gotten = "val" - difference absolute = "diff_abs" relative = "diff_rel; if (diff_rel>0.0001){print "  failed"; exit 1;}else{print "  good"; exit 0;} }' ex=$MU val=$VAL
    if [[ $? -ne 0 ]]; then PASSED=1; fi
  fi
  # overall pass
  if [[ $PASSED -ne 0 ]]; then
    echo "testing kernel values: failed"; exit 1;
  else
    echo "testing kernel values: all good"
  fi
}

# test example
cd $dir

# default setup
# limit time steps for testing
sed -i "s:^NSTEP .*:NSTEP    = 200:" DATA/Par_file
# shortens output interval to avoid timeouts
sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file

# limit time steps for specific examples
# simple mesh example
if [ "$TESTDIR" == "EXAMPLES/applications/meshfem3D_examples/simple_model/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 800:" DATA/Par_file
fi
# tpv5 example
if [ "$TESTDIR" == "EXAMPLES/applications/fault_examples/tpv5/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 1500:" DATA/Par_file
fi
# layered halfspace example
if [ "$TESTDIR" == "EXAMPLES/applications/layered_halfspace/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 200:" DATA/Par_file
fi
# small adjoint example
if [ "$TESTDIR" == "EXAMPLES/applications/small_adjoint_multiple_sources/" ]; then
  # full length
  sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
fi
# socal examples
if [ "$TESTDIR" == "EXAMPLES/applications/meshfem3D_examples/socal1D/" ]; then
  # full length
  sed -i "s:^NSTEP .*:NSTEP    = 840:" DATA/Par_file
  # model setup
  if [ "$TESTID" == "1" ]; then
    # 1D_socal
    sed -i "s:^MODEL .*:MODEL    = 1d_socal:" DATA/Par_file
    rm -f REF_SEIS; ln -s REF_SEIS.1d_socal REF_SEIS
  elif [ "$TESTID" == "2" ]; then
    # 1D_prem
    sed -i "s:^MODEL .*:MODEL    = 1d_prem:" DATA/Par_file
    rm -f REF_SEIS; ln -s REF_SEIS.1d_prem REF_SEIS
  elif [ "$TESTID" == "3" ]; then
    # 1D_cascadia
    sed -i "s:^MODEL .*:MODEL    = 1d_cascadia:" DATA/Par_file
    rm -f REF_SEIS; ln -s REF_SEIS.1d_cascadia REF_SEIS
  else
    # default
    # just continue
    :
  fi
fi
# coupling FK
if [ "$TESTDIR" == "EXAMPLES/applications/small_example_coupling_FK_specfem/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
fi
# elastic halfspace, no absorbing
if [ "$TESTDIR" == "EXAMPLES/applications/homogeneous_halfspace_HEX8_elastic_no_absorbing/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 600:" DATA/Par_file
fi
# waterlayered_halfspace example
if [ "$TESTDIR" == "EXAMPLES/applications/waterlayered_halfspace/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 600:" DATA/Par_file
fi
# tomographic model
if [ "$TESTDIR" == "EXAMPLES/applications/tomographic_model/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
fi
# cavity example
if [ "$TESTDIR" == "EXAMPLES/applications/meshfem3D_examples/cavity/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 2000:" DATA/Par_file
fi
# SEP example
if [ "$TESTDIR" == "EXAMPLES/applications/meshfem3D_examples/sep_bathymetry/" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
fi

# hdf5 i/o example
if [[ "${TEST}" == *"with-hdf5"* ]]; then
  echo
  echo "test run: ${TEST}"
  echo
  sed -i "s:^HDF5_ENABLED .*:HDF5_ENABLED    = .true.:" DATA/Par_file
  sed -i "s:^HDF5_FOR_MOVIES .*:HDF5_FOR_MOVIES    = .true.:" DATA/Par_file
  sed -i "s:^HDF5_IO_NODES .*:HDF5_IO_NODES    = 1:" DATA/Par_file
  # replaces run script
  cp -v run_this_example_HDF5_IO_server.sh run_this_example.sh
fi

# adios
if [ "${ADIOS2}" == "true" ]; then
  # turns on ADIOS
  sed -i "s:^ADIOS_ENABLED .*:ADIOS_ENABLED = .true.:" DATA/Par_file
fi

# GPU
if [ "${GPU}" == "true" ]; then
  # turns on GPU
  sed -i "s:^GPU_ENABLED .*:GPU_ENABLED = .true.:" DATA/Par_file
fi

# use kernel script
if [ "${RUN_KERNEL}" == "true" ]; then
  # use kernel script
  ./run_this_example_kernel.sh
else
  # default script
  ./run_this_example.sh
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# simulation done
echo
echo "simulation done: `pwd`"
echo `date`
echo

# seismogram comparison
if [ "${DEBUG}" == "true" ] || [ "${RUN_KERNEL}" == "true" ]; then
  # no comparisons
  :     # do nothing
else
  my_test
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# kernel test
if [ "${RUN_KERNEL}" == "true" ]; then
  echo
  echo "kernel test directory: `pwd`"
  echo
  my_kernel_test

  # re-run kernel test w/ UNDO_ATT
  # clean up
  rm -rf OUTPUT_FILES/

  # w/ undoatt iteration
  # turns on UNDO_ATTENUATION_AND_OR_PML
  sed -i "s:^UNDO_ATTENUATION_AND_OR_PML .*:UNDO_ATTENUATION_AND_OR_PML = .true.:" DATA/Par_file
  echo
  echo "run kernel w/ UNDO_ATTENUATION_AND_OR_PML"
  echo

  # use kernel script
  ./run_this_example_kernel.sh
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi

  # simulation done
  echo
  echo "simulation done: `pwd`"
  echo `date`
  echo

  # kernel test
  my_kernel_test
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
fi

# cleanup
rm -rf OUTPUT_FILES/
if [ -e DATABASES_MPI ]; then rm -rf DATABASES_MPI/; fi

echo
echo "all good"
echo `date`
echo
