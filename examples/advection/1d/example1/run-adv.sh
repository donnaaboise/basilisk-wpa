#!/usr/bin/env bash

make clean
rm -rf adv
CFLAGS="" make adv.tst
