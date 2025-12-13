//
// Created by Michal Ciesla on 20.10.2022.
//

#ifndef RSA3D_DISCRETEORIENTATIONSSHAPE2_1_H
#define RSA3D_DISCRETEORIENTATIONSSHAPE2_1_H

#include <memory>
#include <vector>

#include "../Shape.h"

/**
 * @brief Class which encapsulated another Shape<2, 1> and samples this particles with gaussian-distribution-drawn
 * orientations.
 *
 * Distribution properies and encapsulated class are set using OrderedShape2_1::initClass. All methods are appropriately
 * delegated to underlying Shape.
 */
class DiscreteOrientationsShape2_1 : public Shape<2, 0>
{
private:
    std::unique_ptr<Shape<2, 1>> underlyingShape;
    static void initAngles(const std::string &args);

public:
	static std::vector<RSAOrientation> allowedOrientations;

	using InitClassFunction = void(*)(const std::string &);

    /**
     * @brief Inits DiscreteOrientationsShape2_1 class, which encapsulated Shape<2, 1> class.
     * @param args Arguments in format:
     *        (n angle1 angle 2 ... angle n)  (encapsulated class parameters)
     * @param baseShapeInitClassFunction pointer to encapsulate Shape<2, 1> initClass function
     */
    static void initClass(const std::string &args, InitClassFunction baseShapeInitClassFunction);
    static Shape<2, 0> *createShape(RND *rnd);

    explicit DiscreteOrientationsShape2_1(std::unique_ptr<Shape<2, 1>> underlyingShape);
    DiscreteOrientationsShape2_1(const DiscreteOrientationsShape2_1 &other);
    DiscreteOrientationsShape2_1 &operator=(const DiscreteOrientationsShape2_1 &other);
    DiscreteOrientationsShape2_1(DiscreteOrientationsShape2_1 &&other) = default;
    DiscreteOrientationsShape2_1 &operator=(DiscreteOrientationsShape2_1 &&other) = default;
    ~DiscreteOrientationsShape2_1() override = default;

    const Vector<2> &getPosition() const override;
    void translate(const Vector<2> &v) override;

    bool overlap(BoundaryConditions<2> *bc, const Shape<2, 0> *s) const override;
    double getVolume(unsigned short dim) const override;
    bool voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition, const Orientation<0> &voxelOrientation,
                     double spatialSize, const Orientation<0> &angularSize) const override;
    double minDistance(const Shape<2, 0> *s) const override;

    std::string toString() const override;
    std::string toPovray() const override;
    std::string toWolfram() const override;

    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;

    Shape<2, 0> *clone() const override;
};


#endif //RSA3D_DISCRETEORIENTATIONSSHAPE2_1_H
