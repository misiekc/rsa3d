//
// Created by PKua on 20.07.18.
//

#ifndef RSA3D_REGULARSOLIDBASE_H
#define RSA3D_REGULARSOLIDBASE_H


#include "../../ConvexShape.h"
#include "../../OverlapStrategyShape.h"
#include "../../../Intersection.h"

class RegularSolidBase : public ConvexShape<3, 0>, public OverlapStrategyShape<3, 0> {
private:
    enum PIResult {
        TRUE,
        FALSE,
        UNKNOWN
    };

    static void printFaceHelperNotebook();
    static void normalizeFacesOrientation();
    static void discoverEdges();
    static void discoverTriangles();
    static void normalizeVolume();
    static void calculateRadia();
    static void discoverAxes();
    static void reportCalculations();
    static intersection::face3D faceFromVertexIdx(const std::vector<size_t> &vertexIdx);
    static void addUniqueAxis(std::vector<Vector<3>> &axes, const Vector<3> &newAxis);
    static void printAxes(const std::vector<Vector<3>> &axes);

    Matrix<3, 3> orientation = Matrix<3, 3>::identity();

    PIResult pointInsideEarlyRejection(const Vector<3> &bcPos) const;
    PIResult pointInsideFace(const Vector<3> &point, const std::vector<Vector<3>> &vertices,
                             const std::vector<Vector<3>> &faceNormals) const;
    bool projectionInsideFace(const Vector<3> &point, const std::vector<Vector<3>> &vertices,
                              const std::vector<size_t> &face, const Vector<3> &faceNormal) const;
    bool pointInsideEdge(const Vector<3> &point, const std::vector<Vector<3>> &vertices) const;
    bool pointInsideVertex(const Vector<3> &point, const std::vector<Vector<3>> &vertices) const;

protected:
    class Edge {
        std::size_t _first{};
        std::size_t _second{};

    public:
        Edge() = default;
        Edge(std::size_t first, std::size_t second) : _first{first}, _second{second} {
            if (_first > _second)   std::swap(_first, _second);
        }
        std::size_t first() const { return _first; }
        std::size_t second() const { return _second; }
        bool operator==(const Edge &other) const { return _first == other._first && _second == other._second; }
    };

    static std::vector<Vector<3>>                   orientedVertices;
    static std::vector<Edge>                        orientedEdges;
    static std::vector<std::vector<std::size_t>>    orientedFaces;
    static std::vector<std::array<std::size_t, 3>>  orientedTriangles;
    static std::vector<Vector<3>>                   orientedFaceNormals;

    static std::vector<Vector<3>> orientedFaceAxes;         // same as face normal, but u,v: u = -v count as one
    static std::vector<Vector<3>> orientedEdgeAxes;         // from vertex to vertex, u = -v as one
    static std::vector<Vector<3>> orientedVertexAxes;       // from center to vertex, u = -v
    static std::vector<Vector<3>> orientedMidedgeAxes;      // from center to midedge, u = -v

    static double normalizeFactor;
    static double circumsphereRadius;
    static double insphereRadius;

    static void initClass(const std::string &attr);

    explicit RegularSolidBase(const Matrix<3, 3> &orientation) : orientation(orientation) {};

    inline std::vector<Vector<3>> applyOrientation(const std::vector<Vector<3>> &vectors) const;
    inline std::vector<Vector<3>> applyPosition(const std::vector<Vector<3>> &vectors) const;
    bool strictPointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation,
                           double orientationRange) const;

public:
    constexpr static double goldRatio = (1 + std::sqrt(5.)) / 2;
    static constexpr double tribonacciConstant = (1 + std::cbrt(19 + 3*std::sqrt(33.)) + std::cbrt(19 - 3*std::sqrt(33.)))/3;

    bool pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation,
                     double orientationRange) const override;
    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
    std::string toWolfram() const override;
    std::vector<std::string> getSupportedStrategies() const override;
    OverlapStrategy<3, 0> *createStrategy(const std::string &name) const override;

    std::vector<Vector<3>> getVertices() const;
    std::vector<Vector<3>> getFaceAxes() const;
    std::vector<Vector<3>> getEdgeAxes() const;
    std::vector<Vector<3>> getVertexAxes() const;
    std::vector<Vector<3>> getMidegdeAxes() const;
    intersection::tri_polyh getTriangles() const;
    intersection::face_polyh getFaces() const;
    std::vector<Vector<3>> getFaceNormals() const;

    const Matrix<3, 3> &getOrientationMatrix() const;
    void setOrientationMatrix(const Matrix<3, 3> &orientation);
    virtual double projectionHalfsize(const Vector<3> &axis) const;
};


#endif //RSA3D_REGULARSOLIDBASE_H
