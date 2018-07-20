//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULARSOLID_H
#define RSA3D_REGULARSOLID_H


#include "../../../Matrix.h"
#include "../../OverlapStrategyShape.h"
#include "../../../Intersection.h"
#include "../../ConvexShape.h"
#include "SATOverlapRS.h"
#include "RegularSolidBase.h"

// CRTP idiom
template <typename SpecificSolid>
class RegularSolid : public RegularSolidBase {
protected:
    static const SATOverlapRS overlapStrategy;   // SATOverlapCB is the default

    explicit RegularSolid(const Matrix<3, 3> &orientation) : RegularSolidBase{orientation} {};

public:
    static void initClass(const std::string &attr);

    bool overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const final;
    Shape<3, 0> *clone() const override;
};

#include "RegularSolid.tpp"

#endif //RSA3D_SOLID_H
