//
// Created by PKua on 03.12.17.
//

#ifndef RSA3D_MINEOVERLAP_H
#define RSA3D_MINEOVERLAP_H

#include "CuboidOverlapStrategy.h"

class MineOverlapCB : public CuboidOverlapStrategy {
private:
    bool checkSegment(const Cuboid *cube, const Vector<3> &point1, const Vector<3> &point2) const;

public:
    bool overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const override;
    std::string getName() const override;

    void runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) const override;
};

#endif //RSA3D_MINEOVERLAP_H
