#!/bin/sh

mpisubmit.bg --stdout "$1"1k.out -n 1 -w 00:01:00 ./"$1" -- 1000 1000
mpisubmit.bg --stdout "$1"2k.out -n 1 -w 00:02:00 ./"$1" -- 2000 2000
mpisubmit.bg --stdout "$1"4k.out -n 1 -w 00:04:00 ./"$1" -- 4000 4000
mpisubmit.bg --stdout "$1"8k.out -n 1 -w 00:08:00 ./"$1" -- 8000 8000
