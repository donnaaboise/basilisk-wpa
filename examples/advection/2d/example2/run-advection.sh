#!/usr/bin/env bash

make clean
rm -rf advection
CC99="mpicc -D_MPI=4" make advection.tst
