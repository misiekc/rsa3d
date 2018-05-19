// (C)PKua 2018
//--------------------------------------------------------------------------------------------

#ifndef _SHAPE_STATIC_SIZES_TEST_H
    #define _SHAPE_STATIC_SIZES_TEST_H
    
#include <ostream> 
#include <vector>
#include "../rsa3d/Parameters.h"
#include "../rsa3d/shapes/OverlapStrategy.h"
#include "utility/ShapePairFactory.h"

/**
 * @brief A simple program testing voxel and neighbour grid cell sizes if they meet the Shape contract.
 *
 * Usage:
 * <blockquote>
 * ./rsa_test shape_sizetest [particle] [attibutes] [max_tries]
 * </blockquote>
 *
 * It simply generates @a max_tries pairs of shapes with random orientation distant by two linear sizes of a voxel or
 * one linear size of neighbour list cell and:
 * <ul>
 * <li> for voxel, checks if a pair always overlaps,
 * <li> for neighbour grid cell, checks if a pair never overlaps.
 * </ul>
 *
 * If those criteria are not met, a conflict occurs and is reported.
 *
 * Next, the program tries a bit smaller voxels and bigger neighbour list cells to see whether the values set are
 * optimal.
 */
namespace shape_sizetest
{
    /**
     * @brief Main function for the test of Shape static sizes.
     *
     * Usage:
     * <blockquote>./rsa_test shape_sizetest [particle] [attibutes] [max_tries]</blockquote>
     */
    int main(int argc, char **argv);
}

#endif  // _SHAPE_STATIC_SIZES_TEST_H
