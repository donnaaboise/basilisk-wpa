#!/usr/bin/env bash

make clean
rm -rf adv
CFLAGS="-events" make adv.tst
