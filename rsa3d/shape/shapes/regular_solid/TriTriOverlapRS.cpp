//
// Created by PKua on 20.07.18.
//

#include "TriTriOverlapRS.h"
#include "RegularSolidBase.h"

bool TriTriOverlapRS::overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const {
    auto &firstSpecific = dynamic_cast<const RegularSolidBase&>(*first);
    auto &secondSpecific = dynamic_cast<const RegularSolidBase&>(*second);

    return intersection::polyh_polyh(firstSpecific.getTriangles(), secondSpecific.getTriangles());
}
