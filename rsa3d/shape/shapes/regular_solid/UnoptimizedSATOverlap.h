//
// Created by PKua on 13.07.18.
//

#ifndef RSA3D_UNOPTIMIZEDSATOVERLAP_H
#define RSA3D_UNOPTIMIZEDSATOVERLAP_H

#include "../../OverlapStrategy.h"

template <typename SpecificSolid>
class UnoptimizedSATOverlap : public OverlapStrategy<3, 0> {
private:
    using interval = std::pair<double, double>;
    using vertices = std::vector<Vector<3>>;

    interval getProjection(const Vector<3> &axis, const vertices &vert) const {
        // Find enpoints of polyhedron projection (multiplied by unknown but const for axis factor)
        interval projInterval = {std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
        for (const auto &v : vert) {
            double proj = v * axis;
            if (proj < projInterval.first)
                projInterval.first = proj;
            if (proj > projInterval.second)
                projInterval.second = proj;
        }
        return projInterval;
    }

    bool isSeparatingAxis(const Vector<3> &axis, const vertices &firstVert, const vertices &secondVert) const {
        interval thisInterval = this->getProjection(axis, firstVert);
        interval otherInterval = this->getProjection(axis, secondVert);

        // TODO epsilon needed
        return std::min(thisInterval.second, otherInterval.second) < std::max(thisInterval.first, otherInterval.first);
    }

public:
    bool overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const override {
        auto firstSpecific = dynamic_cast<const SpecificSolid &>(*first);
        auto secondSpecific = dynamic_cast<const SpecificSolid &>(*second);
        auto firstVertices = firstSpecific.getVertices();
        auto secondVertices = secondSpecific.getVertices();

        // Face axes for this
        for (const auto &face : firstSpecific.getFaceAxes())
            if (this->isSeparatingAxis(face, firstVertices, secondVertices))
                return false;
        // Face axes for other
        for (const auto &face : secondSpecific.getFaceAxes())
            if (this->isSeparatingAxis(face, firstVertices, secondVertices))
                return false;
        // Egde axes cross products
        for (const auto &thisEdge : firstSpecific.getEdgeAxes())
            for (const auto &otherEdge : secondSpecific.getEdgeAxes())
                if (this->isSeparatingAxis(thisEdge ^ otherEdge, firstVertices, secondVertices))
                    return false;

        return true;
    }
};

#endif //RSA3D_UNOPTIMIZEDSATOVERLAP_H
