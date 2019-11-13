//
// Created by PKua on 16.06.18.
//

/**
 * @file ShapeGenerators.h
 * @brief Functions creating Shapes of given position and orientation using ShapeFactory.
 */

#ifndef RSA3D_SHAPEGENERATORS_H
#define RSA3D_SHAPEGENERATORS_H

#include "../../rsa3d/shape/ShapeFactory.h"
#include "../../rsa3d/geometry/Vector.h"
#include "../../rsa3d/utils/Utils.h"
#include <memory>

/**
 * @brief Generates a shape using ShapeFactory of given position and orientation.
 *
 * If provided dimensions are different that ones from ShapeFactory, an exception will be thrown. It can be useful when
 * utilized in tests tailored for specific dimensions (like as2d_extest) which are anyway included in all test
 * executables, even if dimensions are incompatible.
 *
 * A number of orientational degrees of freedom is determined by @a AD and can vary from the real one - note, that some
 * shapes do not support oriented voxels and they pick their orientation by themselves and do not expose it; in that
 * case, the orientation cannot be specified here and will be random.
 *
 * If @a rnd is known not to be used, nullptr can be passed and is a default value of this parameter.
 *
 * @tparam SD spatial dimension of the generated shape
 * @tparam AD angular dimension of the generated shape
 * @param pos position of the generated shape
 * @param orientation orientation of the generated shape
 * @param rnd RND passed to Shape's provided create shape function
 * @return a shape of given position and orientation
 */
template<std::size_t SD = RSA_SPATIAL_DIMENSION, std::size_t AD = RSA_ANGULAR_DIMENSION>
std::unique_ptr<Shape<SD, AD>> generate_shape(const Vector<SD> &pos, const std::array<double, AD> &orientation,
                                              RND *rnd = nullptr) {
    auto shape = dynamic_cast<Shape<SD, AD>*>(ShapeFactory::createShape(rnd));
    if (shape == nullptr) {
        std::ostringstream out;
        out << "Only Shape<" << SD << ", " << AD << "> is supported";
        throw std::runtime_error(out.str());
    }

    shape->translate(pos);
    shape->rotate(orientation);
    return std::unique_ptr<Shape<SD, AD>>(shape);
}

/**
 * @brief Generates a shape using ShapeFactory of given position and random orientation.
 *
 * Note, that some shapes do not support oriented voxels and they pick their orientation randomly by themselves and do
 * not expose it, so orientation will be implicitly random. However, if they do support oriented voxels, the exposed
 * orientational degrees of freedom will be set random explicitly in this method - eventually all real orientational
 * degrees of freedom be random.
 *
 * @tparam SD spatial dimension of the generated shape
 * @tparam AD angular dimension of the generated shape
 * @param pos position of the generated shape
 * @param rnd RND passed to Shape's provided create shape function
 * @return a shape of given position and orientation
 */
template<std::size_t SD = RSA_SPATIAL_DIMENSION, std::size_t AD = RSA_ANGULAR_DIMENSION>
std::unique_ptr<Shape<SD, AD>> generate_randomly_oriented_shape(const Vector<SD> &pos, RND *rnd) {
    double angularSize = Shape<SD, AD>::getAngularVoxelSize();
    std::array<double, AD> orientation{};
    std::for_each(orientation.begin(), orientation.end(), [rnd, angularSize](double &elem) {
        elem = rnd->nextValue() * angularSize;
    });

    return generate_shape<SD, AD>(pos, orientation, rnd);
}


#endif //RSA3D_SHAPEGENERATORS_H