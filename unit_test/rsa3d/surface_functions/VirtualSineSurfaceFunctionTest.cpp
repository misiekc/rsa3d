//
// Created by Piotr Kubala on 25/10/2020.
//

#include <catch.hpp>
#include <iostream>

#include "../../rsa3d/surface_functions/VirtualSineSurfaceFunction.h"

namespace {
    double function_value(const VirtualSineSurfaceFunction &sf, double x) {
        RSAVector v = {{x, 0}};
        sf.fillInLastCoordinate(v);
        return v[1];
    }

    auto function_range(const VirtualSineSurfaceFunction &sf, double x, double dx) {
        RSAVector v = {{x, 0}};
        return sf.calculateValueRange(v, dx);
    }
}

// This test is for visual inspection

/*#include "../../rsa3d/RND.h"

TEST_CASE("VirtualSineSurfaceFunctionTest") {
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

TEST_CASE("VirtualSineSurfaceFunctionTest: bug fixes") {
    // This is for the bug where the argument for maximal value of sine in VirtualSineSurfaceFunction::calculateMinValue
    // was pi instead of pi/2 which yielded improper minimum of virtual sine function
    SECTION("value range bug") {
        VirtualSineSurfaceFunction sf(2, 3.141592653589794, 0.6203504908994);
        double spatialSize = 0.015625;
        RSAVector v{{5.484375, 0}};

        auto[min, max] = sf.calculateValueRange(v, spatialSize);

        CHECK(min == Approx(0.7907537668516279));
        CHECK(max == Approx(0.8821067963959804));
    }
}

// This test is a regression test - the values are taken by using the tested class after visual inspection using
// the above, commented out code
TEST_CASE("VirtualSineSurfaceFunctionTest: values and ranges", "[regression]") {
    SECTION("smooth curve") {
        VirtualSineSurfaceFunction sf(2, 0.2*M_PI, 0.6);

        CHECK(function_value(sf, 2) == Approx(2.532829592264786));
        CHECK(function_value(sf, 6) == Approx(-0.2804049825240839));
        CHECK(function_value(sf, 40) == Approx(0.9552740038267366));

        CHECK(function_range(sf, 3, 0.5).min == Approx(2.329545197473283));
        CHECK(function_range(sf, 3, 0.5).max == Approx(2.532829592264786));
        CHECK(function_range(sf, 2, 1).min == Approx(2.532829592264786));
        CHECK(function_range(sf, 2, 1).max == Approx(2.6));
        CHECK(function_range(sf, 6, 2).min == Approx(-1.4));
        CHECK(function_range(sf, 6, 2).max == Approx(-0.2804049825240839));
        CHECK(function_range(sf, 11, 2).min == Approx(1.989215811185905));
        CHECK(function_range(sf, 11, 2).max == Approx(2.6));
        CHECK(function_range(sf, 1, 8).min == Approx(-1.4));
        CHECK(function_range(sf, 1, 8).max == Approx(2.6));
    }

    SECTION("curve with a little cusp") {
        VirtualSineSurfaceFunction sf(2, 0.4*M_PI, 0.6);

        CHECK(function_value(sf, 2) == Approx(2.259037976638977));
        CHECK(function_value(sf, 3) == Approx(0.3661848922260927));
        CHECK(function_value(sf, 3.74) == Approx(-1.271469875809806));

        CHECK(function_range(sf, 2, 1).min == Approx(0.3661848922260927));
        CHECK(function_range(sf, 2, 1).max == Approx(2.259037976638977));
        CHECK(function_range(sf, 1, 1).min == Approx(2.259037976638977));
        CHECK(function_range(sf, 1, 1).max == Approx(2.6));
        CHECK(function_range(sf, 3, 1).min == Approx(-1.286307093544215));
        CHECK(function_range(sf, 3, 1).max == Approx(0.3661848922260927));
        CHECK(function_range(sf, 9, 1.5).min == Approx(-0.8312116656206461));
        CHECK(function_range(sf, 9, 1.5).max == Approx(2.259037976638976));
        CHECK(function_range(sf, 1, 8).min == Approx(-1.286307093544215));
        CHECK(function_range(sf, 1, 8).max == Approx(2.6));
    }

    SECTION("curve with a big cusp") {
        VirtualSineSurfaceFunction sf(2, M_PI, 0.6);

        CHECK(function_value(sf, 0.75) == Approx(2.55022914280166));
        CHECK(function_value(sf, 1.4) == Approx(1.233334480266489));
        CHECK(function_value(sf, 6.4) == Approx(2.592273266152007));

        CHECK(function_range(sf, 0.7, 0.3).min == Approx(2.372692795266106));
        CHECK(function_range(sf, 0.7, 0.3).max == Approx(2.568565073571581));
        CHECK(function_range(sf, 0.2, 0.7).min == Approx(2.464260158899666));
        CHECK(function_range(sf, 0.2, 0.7).max == Approx(2.6));
        CHECK(function_range(sf, 1.1, 0.7).min == Approx(0.6676929963137455));
        CHECK(function_range(sf, 1.1, 0.7).max == Approx(2.236973863823151));
        CHECK(function_range(sf, 5.9, 0.2).min == Approx(2.236973863823151));
        CHECK(function_range(sf, 5.9, 0.2).max == Approx(2.464260158899665));
        CHECK(function_range(sf, 3.6, 5).min == Approx(0.6676929963137455));
        CHECK(function_range(sf, 3.6, 5).max == Approx(2.6));
    }
}