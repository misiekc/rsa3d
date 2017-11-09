//--------------------------------------------------------------------------------------------
// Test of cuboid intersection algorithm. It uses standard triangle-triangle intersection
// algorithm to check results given by Cuboid::overlap
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _INTTEST_H
    #define _INTTEST_H

namespace cube_inttest
{
    void perform(double sizex, double sizey, double sizez, double box_halfsize, int max_tries);
    void TriTriInt_selftest_run();
}

#endif  // _INTTEST_H
