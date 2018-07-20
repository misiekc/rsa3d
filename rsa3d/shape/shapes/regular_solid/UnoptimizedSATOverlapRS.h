//
// Created by PKua on 13.07.18.
//

#ifndef RSA3D_UNOPTIMIZEDSATOVERLAP_H
#define RSA3D_UNOPTIMIZEDSATOVERLAP_H

#include "../../OverlapStrategy.h"
#include "RegularSolidBase.h"

class UnoptimizedSATOverlapRS : public OverlapStrategy<3, 0> {
private:
    using interval = std::pair<double, double>;
    using vertices = std::vector<Vector<3>>;

    interval getProjection(const Vector<3> &axis, const vertices &vert) const;
    bool isSeparatingAxis(const Vector<3> &axis, const vertices &firstVert, const vertices &secondVert) const;

public:
    bool overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const override;
};

#endif //RSA3D_UNOPTIMIZEDSATOVERLAP_H
