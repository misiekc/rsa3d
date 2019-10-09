//
// Created by PKua on 20.07.18.
//

#ifndef RSA3D_REGULARSOLIDBASE_H
#define RSA3D_REGULARSOLIDBASE_H

#include <memory>
#include <utility>

#include "../../ConvexShape.h"
#include "../../OverlapStrategyShape.h"
#include "../../../geometry/Geometry.h"
#include "../../OrderCalculable.h"

/**
 * @brief Class of common base functionality of concrete regular solid subclasses.
 *
 * The class contains everything which does not require using templates (Shape::clone, for example, requires).
 * Template-related things and also initClass required by ShapeFactory are in RegularSolid subclass.
 */
class RegularSolidBase : public ConvexShape<3, 0>, public OverlapStrategyShape<3, 0>, public OrderCalculable {
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

    /**
     * @brief Helper struct containing all data for the specific RegularSolid.
     */
    struct ShapeData {
        std::vector<Vector<3>> orientedVertices;
        std::vector<Edge> orientedEdges;
        std::vector<std::vector<std::size_t>> orientedFaces;
        std::vector<std::array<std::size_t, 3>> orientedTriangles;
        std::vector<Vector<3>> orientedFaceNormals;

        std::vector<Vector<3>> orientedFaceAxes;         // same as face normal, but u,v: u = -v count as one
        std::vector<Vector<3>> orientedEdgeAxes;         // from vertex to vertex, u = -v as one
        std::vector<Vector<3>> orientedVertexAxes;       // from center to vertex, u = -v
        std::vector<Vector<3>> orientedMidedgeAxes;      // from center to midedge, u = -v

        double normalizeFactor;

        Face3D faceFromVertexIdx(const std::vector<size_t> &vertexIdx) const;
        void printInfo(std::ostream &out) const;
        void printAxes(const std::vector<Vector<3>> &axes, std::ostream &out) const;
    };

private:
    enum PIResult {
        TRUE,
        FALSE,
        UNKNOWN
    };

    static void printFaceHelperNotebook(const ShapeData &shapeData);
    static void normalizeFacesOrientation(ShapeData &shapeDataToComplement);
    static void discoverEdges(ShapeData &shapeDataToComplement);
    static void discoverTriangles(ShapeData &shapeDataToComplement);
    static void normalizeVolume(ShapeData &shapeDataToComplement);
    static void calculateRadia(const ShapeData &shapeData, double &circumsphereRadius, double &insphereRadius);
    static void discoverAxes(ShapeData &shapeDataToComplement);
    static void addUniqueAxis(std::vector<Vector<3>> &axes, const Vector<3> &newAxis);

    Matrix<3, 3> orientation = Matrix<3, 3>::identity();

    PIResult pointInsideFace(const Vector<3> &point, const std::vector<Vector<3>> &vertices) const;
    bool projectionInsideFace(const Vector<3> &point, const std::vector<Vector<3>> &vertices,
                              const std::vector<size_t> &face, const Vector<3> &faceNormal) const;
    bool pointInsideEdge(const Vector<3> &point, const std::vector<Vector<3>> &vertices) const;
    bool pointInsideVertex(const Vector<3> &point, const std::vector<Vector<3>> &vertices) const;

protected:
     /**
      * @brief Static data of the specific RegularSolid shared by all instances of the same RegularSolid subclass.
      *
      * It is not RegularSolid static field, because we want to have coexisting multiple subclasses of RegularSolid to
      * be able to compute order parameters. The ShapeData are passed in the constructor and are prepared earlier by
      * the SpecificClass::calculateStatic and also RegularSolidBase::complementShapeData.
      */
    std::shared_ptr<ShapeData> shapeData;

    /**
     * @brief Using vertices and faces from @a shapeDataToComplement, calculates the rest of ShapeData.
     *
     * @param shapeDataToComplement data containing only vertices and faces, which are to be complemented
     */
    static void complementShapeData(ShapeData &shapeDataToComplement);

    explicit RegularSolidBase(const Matrix<3, 3> &orientation, std::shared_ptr<ShapeData> shapeData)
        : orientation(orientation), shapeData(std::move(shapeData))
     { };

    inline std::vector<Vector<3>> applyOrientation(const std::vector<Vector<3>> &vectors) const;
    inline std::vector<Vector<3>> applyPosition(const std::vector<Vector<3>> &vectors) const;

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
    PolyhedronTriangulation getTriangulation() const;
    PolyhedronFaces getFaces() const;
    std::vector<Vector<3>> getFaceNormals() const;

    const Matrix<3, 3> &getOrientationMatrix() const;
    void setOrientationMatrix(const Matrix<3, 3> &orientation);
    virtual double projectionHalfsize(const Vector<3> &axis) const;
};


#endif //RSA3D_REGULARSOLIDBASE_H
