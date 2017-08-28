//--------------------------------------------------------------------------------------------
// Test of cuboid intersection algorithm. It uses standard triangle-triangle intersection
// algorithm to check results given by Cuboid::overlap
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _INTTEST_H
    #define _INTTEST_H

void CuboidIntTest_run();
void CuboidIntSATTest_run();
void CuboidIntTimeTest_run();
void TriTriInt_selftest_run();

#endif  // _INTTEST_H