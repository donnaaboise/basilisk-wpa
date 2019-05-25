#!/usr/bin/env bash

make clean
rm -rf advection
CFLAGS="" make advection.tst
