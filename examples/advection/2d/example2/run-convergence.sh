#!/usr/bin/env bash

make clean
rm -rf convergence
CFLAGS="" make convergence.tst
