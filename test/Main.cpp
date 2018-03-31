//
// Created by PKua on 24.03.18.
//

#include <cstdlib>
#include "../rsa3d/Utils.h"
#include "CuboidSpeedTest.h"
#include "CuboidIntTest.h"
#include "CuboidPointInsideTest.h"
#include "AnisotropicShape2DExclusionTest.h"
#include "VectorSpeedTest.h"
#include "AnisotropicShape2DExclusionDrawer.h"


int main(int argc, char **argv) {
    if (argc < 2)
        die("No mode param. Aborting.");

    std::string mode(argv[1]);
    if (mode == "cube_speedtest")
        return cube_speedtest::main(argc, argv);
    else if (mode == "cube_inttest")
        return cube_inttest::main(argc, argv);
    else if (mode == "cube_pitest")
        return cube_pitest::main(argc, argv);
    else if (mode == "as2d_extest")
        return as2d_extest::main(argc, argv);
    else if (mode == "vec_speedtest")
        return vec_speedtest::main(argc, argv);
    else if (mode == "as2d_exdrawer")
        return as2d_exdrawer::main(argc, argv);

    std::cerr << "Unknown test mode: " << mode << std::endl;
    return EXIT_FAILURE;
}