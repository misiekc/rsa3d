//--------------------------------------------------------------------------------------------
// Factory generating cuboids from box of given size
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _BOX_FACTORY_H
    #define _BOX_FACTORY_H

#include "ShapePairFactory.h"


class BoxFactory : public ShapePairFactory
{
private:
    double halfsizeX = 0.5;
    double halfsizeY = 0.5;
    double halfsizeZ = 0.5;
    unsigned int no = 0;
    RND rnd;
  
    static BoxFactory * instance;
    
    BoxFactory();
    Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> * randomShape();

public:
    static BoxFactory * getInstance();
    void setBoxSize(double _halfsize_x, double _halfsize_y, double _halfsize_z);
    ShapePair generate();
    std::string getDescription();
};

#endif // _BOX_FACTORY_H
