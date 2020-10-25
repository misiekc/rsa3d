//
// Created by Piotr Kubala on 25/10/2020.
//

#include <catch.hpp>
#include <iostream>

#include "../../rsa3d/surface_functions/VirtualSineSurfaceFunction.h"
#include "../../rsa3d/RND.h"

// This test is for visual inspection

/*TEST_CASE("SineSurfaceFunctionTest") {
    double A = 3;
    double k = 0.7;
    double r = 1.5;
    double max = 20;
    double step = 0.1;
    std::size_t numRectangles = 20;
    VirtualSineSurfaceFunction sf(A, k, r);

    std::cout << "line = Line[{" << std::endl;
    for (double x = 0; x < max; x += step) {
        Vector<2> vector;
        vector[0] = x;
        sf.fillInLastCoordinate(vector);
        std::cout << vector;
        if (x + step < max)
            std::cout << "," << std::endl;
    }
    std::cout << "}];" << std::endl << std::endl;

    RND rnd(12344);
    std::cout << "r = {" << std::endl;
    for (std::size_t i{}; i < numRectangles; i++) {
        double width = rnd.nextValue()*5;
        double x0 = rnd.nextValue() * (max - width);
        auto minMax = sf.calculateValueRange(Vector<2>{{x0, 0}}, width);
        std::cout << "EdgeForm[RandomColor[]],";
        std::cout << "Rectangle[{" << x0 << ", " << minMax.min << "}, {" << (x0 + width) << "," << minMax.max << "}]";
        if (i + 1 < numRectangles)
            std::cout << "," << std::endl;
    }
    std::cout << "};" << std::endl << std::endl;

    std::cout << "Show[Graphics[{line,FaceForm[Transparent],r},Background->White], ";
    std::cout << "Plot[" << A << "*Sin[" << k << "*x],{x,0," << max << "}]]" << std::endl;
}*/