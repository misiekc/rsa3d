// (C)PKua 2018
//--------------------------------------------------------------------------------------------

#ifndef _SHAPE_STATIC_SIZES_TEST_H
    #define _SHAPE_STATIC_SIZES_TEST_H
    
#include <ostream> 
#include <vector>

#include "../rsa3d/Parameters.h"
#include "../rsa3d/shape/OverlapStrategy.h"
#include "utils/ShapePairFactory.h"

/**
 * @brief A simple program testing voxel and neighbour grid cell sizes if they meet the Shape contract.
 *
 * Usage:
 * <blockquote>
 * ./rsa_test shape_sizetest [particle] [attibutes] [max_tries]
 * </blockquote>
 *
 * <strong>It can also be used, assuming correct sizes, to check overlap method correctness. No errors and warning
 * sould then be displayed.</strong>
 *
 * It simply generates @a max_tries pairs of shapes with random orientation distant by two linear sizes of a voxel or
 * one linear size of neighbour list cell and:
 * <ul>
 * <li> for voxel, checks if a pair always overlaps,
 * <li> for neighbour grid cell, checks if a pair never overlaps.
 * </ul>
 * If those criteria are not met, a conflict occurs and is reported.
 *
 * Next, the program tries a bit smaller voxels and bigger neighbour list cells to see whether the values set are
 * optimal. So, then:
 * <ul>
 * <li> for voxel, checks if a pair sometimes is disjunctive,
 * <li> for neighbour grid cell, checks if a pair sometimes overlaps.
 * </ul>
 * If those criteria are not met, either sizes can be altered to be more optimal or too few pairs have been checked.
 * Anyway, a conflict occurs and is reported.
 *
 * On each conflict <strong>the program returns failure exit code 1.</strong> On success the program returns 0.
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
