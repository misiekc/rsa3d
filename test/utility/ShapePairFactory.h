//--------------------------------------------------------------------------------------------
// Interface of factory producing pairs of random cuboids
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _CUBOID_PAIR_FACTORY_H
    #define _CUBOID_PAIR_FACTORY_H

#include <string>
#include <memory>

#include "../../rsa3d/shape/Shape.h"
#include "../../rsa3d/Parameters.h"


/**
 *
 * @brief A class for creating pairs of random Shapes from a specific uniform space.
 *
 * Derived instances can be useful when performing tests on pairs of shapes.
 *
 * @tparam SD spatial dimension of produced shapes
 * @tparam AD angular dimension of produced shapes
 */
template<unsigned short SD, unsigned short AD>
class ShapePairFactory
{
    using shape = Shape<SD, AD>;

protected:
    RND rnd{};

public:
    /**
     * @brief Pair of two shapes. Takes ownership of given shapes and performs deallocation.
     */
    class ShapePair
    {
    private:
        std::shared_ptr<shape> _first;
        std::shared_ptr<shape> _second;

    public:
        ShapePair() = default;
        ShapePair(shape *first, shape *second) : _first(first), _second(second) { }
        ShapePair(std::unique_ptr<shape> first, std::unique_ptr<shape> second) : _first(std::move(first)),
                                                                                 _second(std::move(second)) { }
        shape *first() { return _first.get(); }
        shape *second() { return _second.get(); }
        const shape *first() const { return _first.get(); }
        const shape *second() const { return _second.get(); }
    };
    
    /**
     * @brief Generates a pair of random shapes from a specific uniform space.
     * @return a pair of random shapes from a specific uniform space
     */
    virtual ShapePair generate() = 0;

    /**
     * @brief Returns a brief description of a factory and its parameters.
     * @return a brief description of a factory and its parameters
     */
    virtual std::string getDescription() const = 0;
};

using RSAShapePairFactory = ShapePairFactory<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>;

#endif // _CUBOID_PAIR_FACTORY_H
