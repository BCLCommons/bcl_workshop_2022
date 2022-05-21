#!/bin/bash

./scripts/generate_threaded_seq_from_resfile.sh 274 resfiles/d4.resfile input/sequences_for_threading/d4_native.seq
./scripts/generate_threaded_seq_from_resfile.sh 274 resfiles/d4_d2like.resfile input/sequences_for_threading/d4_d2like.seq
./scripts/generate_threaded_seq_from_resfile.sh 274 resfiles/d4_d3like.resfile input/sequences_for_threading/d4_d3like.seq
./scripts/generate_threaded_seq_from_resfile.sh 274 resfiles/d4_d5like.resfile input/sequences_for_threading/d4_d5like.seq
