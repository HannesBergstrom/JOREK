#!/bin/sh
# Standalone unit test of models/mod_re_kinetic_equilibrium.f90 against stub
# modules (no MPI/HDF5/solver libraries needed). Exercises re_eq_init,
# re_eq_update_labels, re_eq_source and re_eq_source_derivs on a synthetic
# polar grid replicating the re_kin_equil_600 100 keV case, with FP traps
# and array bounds checks enabled.
#
# This directory is deliberately NOT part of the JOREK build (util/ is not
# in the Makefile DIRS): the stubs shadow real module names.
#
# Usage: ./run_test.sh   (from this directory)
set -e
FC=${FC:-gfortran}
# macOS: help gfortran find the system libraries
if [ "$(uname)" = "Darwin" ] && command -v xcrun >/dev/null 2>&1; then
  export LIBRARY_PATH="${LIBRARY_PATH:+$LIBRARY_PATH:}$(xcrun --show-sdk-path)/usr/lib"
fi
# mod_re_kinetic_equilibrium calls LAPACK dgesv (the 'operator' transplant
# variant); link the platform BLAS/LAPACK.
if [ "$(uname)" = "Darwin" ]; then
  LAPACK_LIBS="-F$(xcrun --show-sdk-path)/System/Library/Frameworks -framework Accelerate"
else
  LAPACK_LIBS="-llapack -lblas"
fi
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
$FC -ffree-line-length-none -fcheck=bounds -ffpe-trap=invalid,zero,overflow -g \
    -J "$TMP" stubs.f90 ../../../models/mod_re_kinetic_equilibrium.f90 \
    test_re_eq_unit.f90 -o "$TMP/test_re_eq" $LAPACK_LIBS
( cd "$TMP" && ./test_re_eq )
