//
// Created by PKua on 03.12.17.
//

#ifndef RSA3D_SATOVERLAP_H
#define RSA3D_SATOVERLAP_H

#include "CuboidOverlapStrategy.h"

class SATOverlap : public CuboidOverlapStrategy {
private:
    typedef std::pair<double, double>   interval;

    bool checkSeparatingAxis(const Vector<3> & _axis, Vector<3> * _vert1, Vector<3> * _vert2) const;
    interval getProjection(const Vector<3> & _axis, Vector<3> * _vert) const;

public:
    int overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const override;
    std::string getName() const override;

    void runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) const override;
};


#endif //RSA3D_SATOVERLAP_H
