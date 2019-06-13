#!/usr/bin/env bash

make clean
rm -rf transport
CFLAGS="" make transport.tst
