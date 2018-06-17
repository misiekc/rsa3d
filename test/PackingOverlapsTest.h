//
// Created by PKua on 17.06.18.
//

#ifndef RSA3D_PACKINGOVERLAPSTEST_H
#define RSA3D_PACKINGOVERLAPSTEST_H

#include "../rsa3d/shape/Shape.h"
#include "../rsa3d/Packing.h"

/**
 * @brief A simple program testing voxel packing whether there are any overlaps.
 *
 * Usage:
 * <blockquote>
 * ./rsa_test pack_ovtest [config file] [packing file]
 * </blockquote>
 *
 * It simply checks all pairs of shapes for contact and prints how many overlaps have benn found.
 */
namespace pack_ovtest {
    /**
     * @brief Checks weather there are any overlaps in a @a packing.
     *
     * Calculation are pereformed in parallel by openmp if available.
     * @param packing Packing to be tested
     * @return true, if any overlaps has been found, false otherwise
     */
    bool test_packing_overlaps(const Packing &packing);

    /**
     * @brief Main function for the test of packing overlaps.
     *
     * Usage:
     * <blockquote>./rsa_test pack_ovtest [config file] [packing file]</blockquote>
     */
    int main(int argc, char **argv);
}

#endif //RSA3D_PACKINGOVERLAPSTEST_H
