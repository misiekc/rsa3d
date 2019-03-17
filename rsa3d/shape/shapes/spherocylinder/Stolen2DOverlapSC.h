//
// Created by pkua on 17.03.19.
//

#ifndef RSA3D_STOLEN2DOVERLAPSC_H
#define RSA3D_STOLEN2DOVERLAPSC_H

#include "../../OverlapStrategy.h"
#include "SpheroCylinder2D.h"

/**
 * @brief OverlapStrategy for Spherocylinder<2> using SpheroCylinder2D::overlap method
 */
class Stolen2DOverlapSC : public OverlapStrategy<2, 0> {
public:
    bool overlap(const Shape<2, 0> *first, const Shape<2, 0> *second) const override;
};


#endif //RSA3D_STOLEN2DOVERLAPSC_H
