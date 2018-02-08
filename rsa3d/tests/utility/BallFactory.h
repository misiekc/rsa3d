//--------------------------------------------------------------------------------------------
// Factory generating cuboids from ball of given radius
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _BALL_FACTORY_H
    #define _BALL_FACTORY_H

#include "ShapePairFactory.h"


class BallFactory : public ShapePairFactory
{
private:
    double radius = 0.5;
    unsigned int no = 0;
    RND rnd;
  
    static BallFactory * instance;
    
    BallFactory();
    Shape * randomShape();

public:
    static BallFactory * getInstance();
    void setRadius(double _radius);
    ShapePair generate();
    std::string getDescription();
};

#endif // _BALL_FACTORY_H
