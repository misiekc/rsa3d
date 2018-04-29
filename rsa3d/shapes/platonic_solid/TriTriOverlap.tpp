//
// Created by PKua on 26.04.18.
//

template <typename SpecificSolid>
int TriTriOverlap<SpecificSolid>::overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const {
    auto firstSpecific = dynamic_cast<const SpecificSolid&>(*first);
    auto secondSpecific = dynamic_cast<const SpecificSolid&>(*second);

    return intersection::polyh_polyh(firstSpecific.getTriangles(), secondSpecific.getTriangles());
}