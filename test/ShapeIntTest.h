//--------------------------------------------------------------------------------------------
// Test of cuboid intersection algorithm. It uses standard triangle-triangle intersection
// algorithm to check results given by Cuboid::overlap
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _INTTEST_H
    #define _INTTEST_H
    
#include <ostream> 
#include <vector>
#include "../rsa3d/Parameters.h"
#include "../rsa3d/shapes/OverlapStrategy.h"
#include "utility/ShapePairFactory.h"

namespace shape_inttest
{
    // Struct representing intersection test results
    //--------------------------------------------------------------------------------------------
    struct Results {
        int tries = 0;
        int missed = 0;
        int intersected = 0;
        int disjunct = 0;
        std::vector<ShapePairFactory::ShapePair> missed_dump;
        
        void free_missed_pairs() {
            // Delete missed pairs memory
            for (auto pair : missed_dump)
                pair.free();
        }
    };

    using OverlapStrategy = OverlapStrategy<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>;
    using OverlapStrategyShape = OverlapStrategyShape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION>;

    Results perform(ShapePairFactory *_factory, const OverlapStrategy *firstStrategy, const OverlapStrategy *secondStrategy,
                        int _max_tries);
    void print_results(Results _results);
    void dump_missed_pairs(Results _results, std::ostream & _ostr);

    int main(int argc, char **argv);
}

#endif  // _INTTEST_H
