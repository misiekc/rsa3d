#!/usr/bin/env bash

cd `dirname $0`

../../../cmake-build-release/stat_test/stat_test.3.3 shape_sizetest Kmer3D 2\ 1 10000000
../../../cmake-build-release/stat_test/stat_test.3.3 shape_pitest Kmer3D 2\ 1 10 10000000