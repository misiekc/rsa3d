//
// Created by PKua on 07.02.18.
//

#include "SpheroCylinder2D.h"

int SpheroCylinder2D::pointInside(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) {
    return 0;
}

int SpheroCylinder2D::pointInside(BoundaryConditions *bc, double *da) {
    return AnisotropicShape2D::pointInside(bc, da);
}

double SpheroCylinder2D::getNeighbourListCellSize() {
    return 0;
}

double SpheroCylinder2D::getVoxelSize() {
    return 0;
}

int SpheroCylinder2D::overlap(BoundaryConditions *bc, Shape *s) {
    return 0;
}

double SpheroCylinder2D::getVolume() {
    return 0;
}

std::string SpheroCylinder2D::toString() {
    return Shape::toString();
}

std::string SpheroCylinder2D::toPovray() const {
    return Shape::toPovray();
}

void SpheroCylinder2D::store(std::ostream &f) const {
    Shape::store(f);
}

void SpheroCylinder2D::restore(std::istream &f) {
    Shape::restore(f);
}
