//
// Created by Michal Ciesla on 8/06/2022.
//

#ifndef RSA3D_ALIGNEDRECTANGLE_H
#define RSA3D_ALIGNEDRECTANGLE_H

#include "../Shape.h"
#include "../../BoundaryConditions.h"
#include "Rectangle.h"
#include "../OrderCalculable.h"
#include "../../utils/Utils.h"
#include <vector>

class AlignedRectangle: public Rectangle{

private:
	static std::vector<double> allowedOrientations;

public:
    static void initClass(const std::string &args);

    AlignedRectangle();
    void setAngle(double d) override;
    bool pointInside(BoundaryConditions<2> *bc, const Vector<2> &da, double angleFrom, double angleTo) const override;
};


#endif //RSA3D_ALIGNEDRECTANGLE_H
