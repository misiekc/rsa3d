//
// Created by pkua on 03.09.2019.
//

#ifndef RSA3D_ORDEREDSHAPE2_1_H
#define RSA3D_ORDEREDSHAPE2_1_H

#include <memory>

#include "../Shape.h"

/**
 * @brief Class which encapsulated another Shape<2, 1> and samples this particles with gaussian-distribution-drawn
 * orientations.
 *
 * Distribution properies and encapsulated class are set using OrderedShape2_1::initClass. All methods are appropriately
 * delegated to underlying Shape.
 */
class OrderedShape2_1 : public Shape<2, 0> {
private:
    class RNDUniformRandomBitGenerator {
    private:
        RND &underlyingRnd;

    public:
        explicit RNDUniformRandomBitGenerator(RND &rnd) : underlyingRnd{rnd} { }

        using result_type = unsigned long;

        result_type min() const { return result_type(0); }
        result_type max() const { return std::numeric_limits<result_type>::max(); }
        result_type operator()() {
            return static_cast<result_type>(underlyingRnd.nextValue() * std::numeric_limits<result_type>::max());
        }
    };

    static std::normal_distribution<double> angleDistribution;
    static double preferredAngle;
    static double angleDistributionSigma;

    std::unique_ptr<Shape<2, 1>> underlyingShape;

public:
    using InitClassFunction = void(*)(const std::string &);

    /**
     * @brief Inits OrderedShape2_1 class, which encapsulated Shape<2, 1> class.
     * @param args Arguments in format:
     *        (angle distribution average) (angle distribution sigma) (encapsulated class parameters)
     * @param baseShapeInitClassFunction pointer to encapsulate Shape<2, 1> initClass function
     */
    static void initClass(const std::string &args, InitClassFunction baseShapeInitClassFunction);
    static Shape<2, 0> *createShape(RND *rnd);

    explicit OrderedShape2_1(std::unique_ptr<Shape<2, 1>> underlyingShape);
    OrderedShape2_1(const OrderedShape2_1 &other);
    OrderedShape2_1 &operator=(const OrderedShape2_1 &other);
    OrderedShape2_1(OrderedShape2_1 &&other) = default;
    OrderedShape2_1 &operator=(OrderedShape2_1 &&other) = default;
    ~OrderedShape2_1() override = default;

    const Vector<2> &getPosition() const override;
    void translate(const Vector<2> &v) override;

    bool overlap(BoundaryConditions<2> *bc, const Shape<2, 0> *s) const override;
    double getVolume(unsigned short dim) const override;
    bool voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition, const Orientation<0> &voxelOrientation,
                     double spatialSize, double angularSize) const override;
    double minDistance(const Shape<2, 0> *s) const override;

    std::string toString() const override;
    std::string toPovray() const override;
    std::string toWolfram() const override;

    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;

    Shape<2, 0> *clone() const override;
};


#endif //RSA3D_ORDEREDSHAPE2_1_H
