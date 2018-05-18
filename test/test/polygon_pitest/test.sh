#!/usr/bin/env bash

cd `dirname $0`

../../../rsa_test as2d_extest test1.txt test.nb
../../../rsa_test as2d_extest test2.txt test.nb
../../../rsa_test as2d_extest test3.txt test.nb
../../../rsa_test as2d_extest test4.txt test.nb
../../../rsa_test as2d_extest test5.txt test.nb
../../../rsa_test as2d_extest test6.txt test.nb
../../../rsa_test as2d_extest test7.txt test.nb
../../../rsa_test as2d_extest test8.txt test.nb
../../../rsa_test as2d_extest test9.txt test.nb
../../../rsa_test as2d_extest test10.txt test.nb
../../../rsa_test as2d_extest test11.txt test.nb
../../../rsa_test as2d_extest test12.txt test.nb
../../../rsa_test as2d_extest test13.txt test.nb
../../../rsa_test as2d_extest test14.txt test.nb

mathematica test.nb &
