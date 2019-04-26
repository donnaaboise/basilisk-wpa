#!/usr/bin/env bash

make clean
rm -rf swe
CFLAGS="" make swe.tst
