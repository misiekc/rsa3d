//
// Created by PKua on 26.04.18.
//

#ifndef RSA3D_TRITRIOVERLAP_H
#define RSA3D_TRITRIOVERLAP_H


#include "../../OverlapStrategy.h"

class TriTriOverlapRS : public OverlapStrategy<3, 0> {
public:
    bool overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const override;
};

#endif //RSA3D_TRITRIOVERLAP_H