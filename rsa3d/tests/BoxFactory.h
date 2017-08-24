//--------------------------------------------------------------------------------------------
// Factory generating cuboids from box of given size
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include "CuboidPairFactory.h"


class BoxFactory : public CuboidPairFactory
{
private:
    double halfsizeX = 0.5;
    double halfsizeY = 0.5;
    double halfsizeZ = 0.5;

public:
    void setSize(double _halfsize_x, double _halfsize_y, double _halfsize_z);
    CuboidPair generate();
};
