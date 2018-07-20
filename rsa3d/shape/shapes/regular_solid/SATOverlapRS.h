//
// Created by PKua on 26.04.18.
//

#ifndef RSA3D_SATOVERLAP_H
#define RSA3D_SATOVERLAP_H


#include "../../OverlapStrategy.h"
#include "RegularSolidBase.h"

class SATOverlapRS : public OverlapStrategy<3, 0> {
private:
    bool isSeparatingAxis(const Vector<3> &axis, const RegularSolidBase &first, const RegularSolidBase &second,
                          const Vector<3> &distance) const;

public:
    bool overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const override;
};


#endif //RSA3D_SATOVERLAP_H
