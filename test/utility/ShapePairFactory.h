//--------------------------------------------------------------------------------------------
// Interface of factory producing pairs of random cuboids
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _CUBOID_PAIR_FACTORY_H
    #define _CUBOID_PAIR_FACTORY_H

#include <string>
#include <memory>

#include "../../rsa3d/Shape.h"
#include "../../rsa3d/Parameters.h"


/**
 * @brief A class for creating pairs of random Shapes from a specific uniform space.
 *
 * Derived instances can be useful when performing tests on pairs of shapes.
 */
class ShapePairFactory
{
protected:
    RND rnd{};

    std::array<double, RSA_ANGULAR_DIMENSION> randomOrientation() {
        std::array<double, RSA_ANGULAR_DIMENSION> orientation{};
        std::for_each(orientation.begin(), orientation.end(), [this](double &elem) {
            elem = this->rnd.nextValue() * RSAShape::getVoxelAngularSize();
        });
        return orientation;
    };

public:
    /**
     * @brief Pair of two shapes. Takes ownership of given shapes and perform deallocation.
     */
    class ShapePair
    {
    private:
        std::shared_ptr<RSAShape> _first;
        std::shared_ptr<RSAShape> _second;

    public:
        ShapePair(RSAShape *first, RSAShape *second) : _first(first), _second(second) { }
        RSAShape *first() { return _first.get(); }
        RSAShape *second() { return _second.get(); }
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

#endif // _CUBOID_PAIR_FACTORY_H
