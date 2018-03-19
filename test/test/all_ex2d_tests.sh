#!/usr/bin/env bash

cd `dirname $0` # change working dir to script dir

./spheroc_inttest/test.sh
./spheroc_pitest/test.sh
./ellipse_inttest/test.sh
./ellipse_pitest/test.sh
./rectangle_inttest/test.sh
./rectangle_pitest/test.sh
