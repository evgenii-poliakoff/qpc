#!/bin/bash

export VERSION=$(cat $SRC_DIR/VERSION)

mkdir -p  $SRC_DIR/build
pushd $SRC_DIR/build

export FC=gfortran
export FFLAGS="-Ofast -ffree-line-length-512"

$PYTHON -m numpy.f2py -c --quiet -m _fastmul $SRC_DIR/src/qpc/fastmul.f90 --backend meson
$PYTHON -m numpy.f2py -c --quiet -m _evolution $SRC_DIR/src/qpc/evolution.f90 --backend meson
$PYTHON -m numpy.f2py -c --quiet -m _evolution_chained2 $SRC_DIR/src/qpc/evolution_chained2.f90 --backend meson
$PYTHON -m numpy.f2py -c --quiet -m _evolution_chained2_kicked $SRC_DIR/src/qpc/evolution_chained2_kicked.f90 --backend meson
$PYTHON -m numpy.f2py -c --quiet -m _fastmuli $SRC_DIR/src/qpc/fastmuli.f90 --backend meson
$PYTHON -m numpy.f2py -c --quiet -m _local_op $SRC_DIR/src/qpc/local_op.f90 --backend meson

cp _fastmul.*.so $SRC_DIR/qpc
cp _evolution.*.so $SRC_DIR/qpc
cp _evolution_chained2.*.so $SRC_DIR/qpc
cp _evolution_chained2_kicked.*.so $SRC_DIR/qpc
cp _fastmuli.*.so $SRC_DIR/qpc
cp _local_op.*.so $SRC_DIR/qpc

popd

$PYTHON setup.py install