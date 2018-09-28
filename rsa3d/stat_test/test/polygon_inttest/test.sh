#!/usr/bin/env bash

cd `dirname $0`

../../../rsa_test as2d_extest test1.txt test.nb
../../../rsa_test as2d_extest test4.txt test.nb
../../../rsa_test as2d_extest test7.txt test.nb

mathematica test.nb &
