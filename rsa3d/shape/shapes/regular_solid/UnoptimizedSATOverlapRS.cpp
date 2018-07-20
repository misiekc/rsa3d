//
// Created by PKua on 20.07.18.
//

#include "UnoptimizedSATOverlapRS.h"

UnoptimizedSATOverlapRS::interval UnoptimizedSATOverlapRS::getProjection(const Vector<3> &axis,
                                                                         const std::vector <Vector<3>> &vert) const {
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

bool UnoptimizedSATOverlapRS::isSeparatingAxis(const Vector<3> &axis,
                                               const UnoptimizedSATOverlapRS::vertices &firstVert,
                                               const UnoptimizedSATOverlapRS::vertices &secondVert) const {
    interval thisInterval = this->getProjection(axis, firstVert);
    interval otherInterval = this->getProjection(axis, secondVert);

    // TODO epsilon needed
    return std::min(thisInterval.second, otherInterval.second) < std::max(thisInterval.first, otherInterval.first);
}

bool UnoptimizedSATOverlapRS::overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const {
    auto &firstSpecific = dynamic_cast<const RegularSolidBase &>(*first);
    auto &secondSpecific = dynamic_cast<const RegularSolidBase &>(*second);
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
