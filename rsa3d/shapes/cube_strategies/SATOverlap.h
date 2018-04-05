//
// Created by PKua on 03.12.17.
//

#ifndef RSA3D_SATOVERLAP_H
#define RSA3D_SATOVERLAP_H

#include "OverlapStrategy.h"

class SATOverlap : public OverlapStrategy {
private:
    typedef std::pair<double, double>   interval;

    bool checkSeparatingAxis(const Vector<3> & _axis, Vector<3> * _vert1, Vector<3> * _vert2) const;
    interval getProjection(const Vector<3> & _axis, Vector<3> * _vert) const;

public:
    bool overlap(const Cuboid *cube1, const Cuboid *cube2, BoundaryConditions *bc) override;
    std::string getName() override;

    void runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) override;
};


#endif //RSA3D_SATOVERLAP_H
