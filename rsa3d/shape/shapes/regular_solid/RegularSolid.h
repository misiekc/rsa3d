//
// Created by PKua on 21.04.18.
//

#ifndef RSA3D_REGULARSOLID_H
#define RSA3D_REGULARSOLID_H


#include "../../../Matrix.h"
#include "../../OverlapStrategyShape.h"
#include "../../../Intersection.h"
#include "../../ConvexShape.h"

// CRTP idiom
template <typename SpecificSolid>
class RegularSolid : public ConvexShape<3, 0>, public OverlapStrategyShape<3, 0> {
protected:
    using interval = std::pair<double, double>;

private:
    Matrix<3, 3> orientation = Matrix<3, 3>::identity();

    static void discoverTriangles();
    static void normalizeVolume();
    static void calculateRadia();
    static void discoverAxes();
    static void addUniqueAxis(std::vector<Vector<3>> &axes, const Vector<3> &newAxis);

    interval getProjection(const Vector<3> & axis) const;

protected:
    static std::vector<Vector<3>> orientedVertices;
    static std::vector<std::vector<std::size_t>> orientedFaces;
    static std::vector<std::array<std::size_t, 3>> orientedTriangles;
    static std::vector<Vector<3>> orientedFaceAxes;
    static std::vector<Vector<3>> orientedEdgeAxes;
    static std::vector<Vector<3>> orientedVertexAxes;
    static std::vector<Vector<3>> orientedMidedgeAxes;

    static double normalizeFactor;
    static double circumsphereRadius;
    static double insphereRadius;

    explicit RegularSolid(const Matrix<3, 3> &orientation) : orientation(orientation) {};

    inline std::vector<Vector<3>> applyOrientation(const std::vector<Vector<3>> &vectors) const;
    inline std::vector<Vector<3>> applyPosition(const std::vector<Vector<3>> &vectors) const;
    void setOrientationMatrix(const Matrix<3, 3> &orientation);

    bool isSeparatingAxisUnoptimized(const Vector<3> &axis, const SpecificSolid &other, const Vector<3> &distance) const;

public:
    static void initClass(const std::string &attr);

    bool overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const final;
    bool pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation,
                     double orientationRange) const override;
    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
    std::string toWolfram() const override;
    Shape<3, 0> *clone() const override;

    std::vector<std::string> getSupportedStrategies() const override;
    OverlapStrategy<3, 0> *createStrategy(const std::string &name) const override;

    std::vector<Vector<3>> getVertices() const;
    std::vector<Vector<3>> getFaceAxes() const;
    std::vector<Vector<3>> getEdgeAxes() const;
    std::vector<Vector<3>> getVertexAxes() const;
    std::vector<Vector<3>> getMidegdeAxes() const;
    intersection::tri_polyh getTriangles() const;
    intersection::face_polyh getFaces() const;

    const Matrix<3, 3> &getOrientationMatrix() const;

    /* double projectionHalfsize(const Vector<3> &axis) const; */   /* CRTP pure virtual */
    bool isSeparatingAxis(const Vector<3> &axis, const SpecificSolid &other,
                          const Vector<3> &distance) const;         /* CRTP virtual */
    static void printFaceHelperNotebook();
};

#include "RegularSolid.tpp"

#endif //RSA3D_SOLID_H
