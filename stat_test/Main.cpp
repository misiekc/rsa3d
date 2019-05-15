//
// Created by PKua on 24.03.18.
//

#include <cstdlib>
#include "../rsa3d/utils/Utils.h"
#include "CuboidSpeedTest.h"
#include "ShapeOverlapTest.h"
#include "ShapePointInsideTest.h"
#include "ShapeStaticSizesTest.h"
#include "AnisotropicShape2DExclusionTest.h"
#include "VectorSpeedTest.h"
#include "AnisotropicShape2DExclusionDrawer.h"
#include "RectangleCase.h"
#include "PackingOverlapsTest.h"
#include "OrderParamTest.h"
#include "ShapeBCTest.h"


int main(int argc, char **argv) {
    if (argc < 2)
        die("No mode param. Aborting.");

    std::string mode(argv[1]);
    if (mode == "cube_speedtest")
        return cube_speedtest::main(argc, argv);
    else if (mode == "shape_ovtest")
        return shape_ovtest::main(argc, argv);
    else if (mode == "shape_pitest")
        return shape_pitest::main(argc, argv);
    else if (mode == "shape_sizetest")
        return shape_sizetest::main(argc, argv);
    else if (mode == "shape_bctest")
        return shape_bctest::main(argc, argv);
    else if (mode == "as2d_extest")
        return as2d_extest::main(argc, argv);
    else if (mode == "vec_speedtest")
        return vec_speedtest::main(argc, argv);
    else if (mode == "as2d_exdrawer")
        return as2d_exdrawer::main(argc, argv);
    else if (mode == "rectangle_case")
    	return rectanglecase::main(argc, argv);
    else if (mode == "pack_ovtest")
        return pack_ovtest::main(argc, argv);
    else if (mode == "order_param_test")
        return order_param_test::main(argc, argv);

    std::cerr << "Unknown test mode: " << mode << std::endl;
    return EXIT_FAILURE;
}
