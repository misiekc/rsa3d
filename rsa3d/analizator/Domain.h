//
// Created by ciesla on 12/27/22.
//

#ifndef RSA3D_DOMAIN_H
#define RSA3D_DOMAIN_H

#include <unordered_set>
#include "../shape/Shape.h"
#include "../NeighbourGrid.h"
#include "../shape/shapes/Rectangle.h"
#include "../shape/shapes/DiscreteOrientationsShape2_1.h"
#include "../Parameters.h"
#include "../boundary_conditions/FreeBC.h"
#include "../boundary_conditions/PeriodicBC.h"

class Domain {

private:
    std::unordered_set<const RSAShape *> shapes;
    NeighbourGrid<const RSAShape> *ng;
    double orientation{};
    double rectangleSides[2] = {Rectangle::longer, Rectangle::shorter};
    double surfaceSize;
    unsigned short int surfaceDimension;


    void addShape(const RSAShape *shape);

public:

    explicit Domain(double surfaceSize, unsigned short int surfaceDimension, const RSAShape *seed);

    ~Domain();

    bool testAndAddShape(const RSAShape *shape, RSABoundaryConditions *bc);

    size_t size() const;
    double getOrientation() const;
    std::unordered_set<const RSAShape *> getShapes() const;
};


#endif //RSA3D_DOMAIN_H
