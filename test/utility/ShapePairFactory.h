//--------------------------------------------------------------------------------------------
// Interface of factory producing pairs of random cuboids
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _CUBOID_PAIR_FACTORY_H
    #define _CUBOID_PAIR_FACTORY_H

#include <string>

#include "../../rsa3d/Shape.h"
#include "../../rsa3d/Parameters.h"


/**
 * @brief A class for creating pairs of random Shapes from a specific uniform space.
 *
 * Derived instances can be useful when performing tests on pairs of shapes.
 */
class ShapePairFactory
{
public:
    using RSAShape = Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>;

    /**
     * @brief Pair of two shapes. Provides also self-cleaning utility.
     */
    struct ShapePair
    {
        RSAShape *first;
        RSAShape *second;

        ShapePair(RSAShape *first, RSAShape *second) : first(first), second(second) { }
        
        /**
         * @brief Deallocates memory for stored shapes.
         */
        inline void free() {
            delete first;
            delete second;
        }
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
