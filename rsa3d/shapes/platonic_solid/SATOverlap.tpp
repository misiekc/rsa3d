//
// Created by PKua on 26.04.18.
//

template <typename SpecificSolid>
int SATOverlap<SpecificSolid>::overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const {
    auto firstSpecific = dynamic_cast<const SpecificSolid&>(*first);
    auto secondSpecific = dynamic_cast<const SpecificSolid&>(*second);
    Vector<3> distance = Vector<3>(secondSpecific.getPosition()) - Vector<3>(firstSpecific.getPosition());

    // Face axes for this
    for (const auto &face : firstSpecific.getFaceAxes())
        if (firstSpecific.isSeparatingAxis(face, secondSpecific, distance))
            return 0;
    // Face axes for other
    for (const auto &face : secondSpecific.getFaceAxes())
        if (firstSpecific.isSeparatingAxis(face, secondSpecific, distance))
            return 0;
    // Egde axes cross products
    for (const auto &thisEdge : firstSpecific.getEdgeAxes())
        for (const auto &otherEdge : secondSpecific.getEdgeAxes())
            if (firstSpecific.isSeparatingAxis(thisEdge ^ otherEdge, secondSpecific, distance))
                return 0;

    return 1;
}
