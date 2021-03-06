//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULARSOLID_H
#define RSA3D_REGULARSOLID_H


#include "../../../geometry/Matrix.h"
#include "../../OverlapStrategyShape.h"
#include "../../../geometry/Geometry.h"
#include "../../ConvexShape.h"
#include "SATOverlapRS.h"
#include "RegularSolidBase.h"


/**
 * @brief Class complementing RegularSolidBase with stuff requiring templates and specific to concrete solids.
 *
 * It uses CRTP idiom - enables concrete solids to provide "polymorphic" calculateStatic function, an overlap strategy
 * of choice and Platonic solid of which symmetry to use to calculate order parameters (see any specific class how it
 * works). It takes care of the rest - registering ShapeStaticInfo, cloning and
 * passing common instance of ShapeData (see RegularSolidBase::shapeData for explanation).
 *
 * @tparam SpecificSolid CRTP-idiom-static-polymorphic class
 */
template <typename SpecificSolid>
class RegularSolid : public RegularSolidBase {
private:
    static std::shared_ptr<ShapeData> staticShapeData;

    static bool isSpecificSolidBorrowingOrderParameters();

protected:
    /**
     * @brief Platonic solid which should be used to calculate order parameters (Tetrahedron, Octahedron or
     * Icosahedron).
     *
     * Default value is only placeholder, most classes (apart from the above three) should "override" it (see any but
     * those three to know how it works).
     */
    using SymmetryPlatonicSolid = SpecificSolid;

    /**
     * @brief Default overlapStrategy. SpecificSolid classes can "override" this field with different strategy class
     * (see for example how Tetrahedron does this).
     */
    static const SATOverlapRS overlapStrategy;   // SATOverlapRS is the default

    explicit RegularSolid(const Matrix<3, 3> &orientation)
        : RegularSolidBase(orientation, RegularSolid::staticShapeData)
    { };

public:
    static void initClass(const std::string &attr);

    bool overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const final;
    Shape<3, 0> *clone() const override;

    std::unique_ptr<RegularSolidBase> createSymmetryShape() const override;

    std::vector<double> calculateOrder(const OrderCalculable *other) const override;
};

#include "RegularSolid.tpp"

#endif //RSA3D_SOLID_H
