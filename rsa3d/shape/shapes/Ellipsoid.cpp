//
// Created by PKua on 17.07.18.
//

#include "Ellipsoid.h"

double Ellipsoid::a{};
double Ellipsoid::b{};
double Ellipsoid::c{};

void Ellipsoid::initClass(const std::string &attr) {
    std::istringstream attrStream(attr);

    attrStream >> a >> b >> c;
    if (!attrStream)    throw std::runtime_error("Wrong attr format");
    // Axes of increasing size a <= b <= c - the largest axis lies on z-axis
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);

    normalizeVolume();

    Shape::setNeighbourListCellSize(2*c);
    Shape::setVoxelSpatialSize(2*a/std::sqrt(3));

    Shape::setCreateShapeImpl([](RND *rnd) -> Shape* {
        return new Ellipsoid(Matrix<3, 3>::rotation(
                2*M_PI*rnd->nextValue(),
                std::asin(2*rnd->nextValue() - 1),
                2*M_PI*rnd->nextValue()));
    });
}

void Ellipsoid::normalizeVolume() {
    double volume = 4./3*M_PI*a*b*c;
    double factor = 1/std::cbrt(volume);
    a *= factor;
    b *= factor;
    c *= factor;
}

bool Ellipsoid::overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const {
    return false;
}

bool Ellipsoid::pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation,
                            double orientationRange) const {
    return false;
}

void Ellipsoid::store(std::ostream &f) const {
    Shape::store(f);
    double d;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
            d = this->orientation(i, j);
            f.write((char *)(&d), sizeof(double));
        }
    }
}

void Ellipsoid::restore(std::istream &f) {
    Shape::restore(f);
    double d;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
            f.read((char *)&d, sizeof(double));
            this->orientation(i, j) = d;
        }
    }
}

Shape<3, 0> *Ellipsoid::clone() const {
    return new Ellipsoid(*this);
}
