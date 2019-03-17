//
// Created by pkua on 17.03.19.
//

#ifndef RSA3D_OWNOVERLAPSC_H
#define RSA3D_OWNOVERLAPSC_H

#include "../../OverlapStrategy.h"
#include "../../../FreeBC.h"

template <unsigned short DIM>
class Spherocylinder;

/**
 * @brief Dummy OverlapStrategy using the very own Spherocylinder::overlap method
 *
 * @tparam DIM Spherocylinder dimension
 */
template <unsigned short DIM>
class OwnOverlapSC : public OverlapStrategy<DIM, 0> {
public:
    bool overlap(const Shape<DIM, 0> *first, const Shape<DIM, 0> *second) const override {
        auto &firstSpecific = dynamic_cast<const Spherocylinder<DIM> &>(*first);
        auto &secondSpecific = dynamic_cast<const Spherocylinder<DIM> &>(*second);

        FreeBC<DIM> bc;
        return firstSpecific.overlap(&bc, &secondSpecific);
    }
};

#endif //RSA3D_OWNOVERLAPSC_H
