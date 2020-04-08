//
// Created by Piotr Kubala on 26/03/2020.
//

#include <catch.hpp>

#include "../../../../rsa3d/shape/shapes/RegularDiskopolygon.h"
#include "../../../../rsa3d/shape/shapes/spherocylinder/SpheroCylinder2D.h"
#include "../../../matchers/VectorApproxMatcher.h"

TEST_CASE("RegularDiskopolygonAttributes") {
    SECTION("standard") {
        RegularDiskopolygonAttributes attr("standard 3 4 5"); // triangle, side 4, radius 5

        REQUIRE(attr.getNSides() == 3);
        REQUIRE(attr.getSideLength() == Approx(0.3316471184014903));
        REQUIRE(attr.getRadius() == Approx(0.4145588980018629));
        REQUIRE(attr.getHeight() == Approx(0.09573827654253206));
        REQUIRE(attr.getHalfDiagonal() == Approx(0.1914765530850641));
    }

    SECTION("short") {
        RegularDiskopolygonAttributes attr("short 3 2"); // triangle, smoothing radius 2, implicit circumradius 1

        REQUIRE(attr.getNSides() == 3);
        REQUIRE(attr.getSideLength() == Approx(0.3516703029486882));
        REQUIRE(attr.getRadius() == Approx(0.4060738881468448));
        REQUIRE(attr.getHeight() == Approx(0.1015184720367112));
        REQUIRE(attr.getHalfDiagonal() == Approx(0.2030369440734224));
    }

    SECTION("validation") {
        SECTION("malformed") {
            REQUIRE_THROWS_WITH(RegularDiskopolygonAttributes("skrewed 1 2 3"), Catch::Contains("Supported"));
            REQUIRE_THROWS_WITH(RegularDiskopolygonAttributes("standard 3 2"), Catch::Contains("Malformed"));
            REQUIRE_THROWS_WITH(RegularDiskopolygonAttributes("standard 3"), Catch::Contains("Malformed"));
            REQUIRE_THROWS_WITH(RegularDiskopolygonAttributes("standard 3 2 a"), Catch::Contains("Malformed"));
            REQUIRE_THROWS_WITH(RegularDiskopolygonAttributes("standard 3 a 1"), Catch::Contains("Malformed"));
            REQUIRE_THROWS_WITH(RegularDiskopolygonAttributes("standard a 2 1"), Catch::Contains("Malformed"));
            REQUIRE_THROWS_WITH(RegularDiskopolygonAttributes("short 3"), Catch::Contains("Malformed"));
            REQUIRE_THROWS_WITH(RegularDiskopolygonAttributes("short 3 a"), Catch::Contains("Malformed"));
            REQUIRE_THROWS_WITH(RegularDiskopolygonAttributes("short a 3"), Catch::Contains("Malformed"));
        }

        SECTION("mathematically incorrect") {
            REQUIRE_THROWS(RegularDiskopolygon::initClass("standard 1 2 3"));
            REQUIRE_THROWS(RegularDiskopolygon::initClass("standard 2 2 3"));
            REQUIRE_THROWS(RegularDiskopolygon::initClass("standard 3 1 0"));
            REQUIRE_NOTHROW(RegularDiskopolygon::initClass("standard 3 0 1"));
            REQUIRE_THROWS(RegularDiskopolygon::initClass("short 1 2"));
            REQUIRE_THROWS(RegularDiskopolygon::initClass("short 2 2"));
            REQUIRE_THROWS(RegularDiskopolygon::initClass("short 3 0"));
        }
    }
}

TEST_CASE("RegularDiskopolygonAttributes: initClass") {
    SECTION("shape info") {
        RegularDiskopolygon::initClass("standard 3 4 5");

        auto shapeInfo = Shape<2, 1>::getShapeStaticInfo();
        REQUIRE(shapeInfo.getCircumsphereRadius() == Approx(0.6060354510869270));
        REQUIRE(shapeInfo.getInsphereRadius() == Approx(0.5102971745443950));
        REQUIRE(shapeInfo.getVoxelAngularSize() == Approx(2.094395102393195));
    }

    SECTION("spherocylinder parameters") {
        RegularDiskopolygon::initClass("standard 3 4 5");

        REQUIRE(SpheroCylinder2D::getRadius() == Approx(0.4145588980018629));
        REQUIRE(SpheroCylinder2D::getDistance() == Approx(0.3316471184014903));
    }
}

TEST_CASE("RectangularBounding") {
    SECTION("single point") {
        RectangularBounding rb;

        rb.addPoint({{1, 2}});

        REQUIRE(rb.getBottomLeft() == Vector<2>{{1, 2}});
        REQUIRE(rb.getTopRight() == Vector<2>{{1, 2}});
    }

    SECTION("two points - all points check") {
        RectangularBounding rb;

        rb.addPoint({{1, 2}});
        rb.addPoint({{3, 4}});

        REQUIRE(rb.getBottomLeft() == Vector<2>{{1, 2}});
        REQUIRE(rb.getBottomRight() == Vector<2>{{3, 2}});
        REQUIRE(rb.getTopLeft() == Vector<2>{{1, 4}});
        REQUIRE(rb.getTopRight() == Vector<2>{{3, 4}});
    }

    SECTION("3 points - 3rd inside") {
        RectangularBounding rb;

        rb.addPoint({{3, 4}});
        rb.addPoint({{1, 2}});
        rb.addPoint({{2, 3}});

        REQUIRE(rb.getBottomLeft() == Vector<2>{{1, 2}});
        REQUIRE(rb.getTopRight() == Vector<2>{{3, 4}});
    }

    SECTION("3 points - 3rd outside") {
        RectangularBounding rb;

        rb.addPoint({{3, 4}});
        rb.addPoint({{1, 2}});

        SECTION("up") {
            rb.addPoint({{2, 5}});

            REQUIRE(rb.getBottomLeft() == Vector<2>{{1, 2}});
            REQUIRE(rb.getTopRight() == Vector<2>{{3, 5}});
        }

        SECTION("down") {
            rb.addPoint({{2, -5}});

            REQUIRE(rb.getBottomLeft() == Vector<2>{{1, -5}});
            REQUIRE(rb.getTopRight() == Vector<2>{{3, 4}});
        }

        SECTION("left") {
            rb.addPoint({{-1, 3}});

            REQUIRE(rb.getBottomLeft() == Vector<2>{{-1, 2}});
            REQUIRE(rb.getTopRight() == Vector<2>{{3, 4}});
        }

        SECTION("right") {
            rb.addPoint({{8, 3}});

            REQUIRE(rb.getBottomLeft() == Vector<2>{{1, 2}});
            REQUIRE(rb.getTopRight() == Vector<2>{{8, 4}});
        }
    }

    SECTION("expansion") {
        RectangularBounding rb;
        rb.addPoint({{1, 2}});
        rb.addPoint({{3, 4}});

        rb.expand(3);

        REQUIRE(rb.getBottomLeft() == Vector<2>{{1, 2}});
        REQUIRE(rb.getTopRight() == Vector<2>{{6, 7}});
    }

    SECTION("translation") {
        RectangularBounding rb;
        rb.addPoint({{1, 2}});
        rb.addPoint({{3, 4}});

        rb.translate(Vector<2>{{-1, 1}});

        REQUIRE(rb.getBottomLeft() == Vector<2>{{0, 3}});
        REQUIRE(rb.getTopRight() == Vector<2>{{2, 5}});
    }
}

TEST_CASE("RectangularBoundingBuilder: building for arch") {
    SECTION("only first quarter") {
        RectangularBounding bounding = RectangularBoundingBuilder::forArch({{0, 1}}, M_PI / 6, M_PI / 4);

        REQUIRE_THAT(bounding.getBottomLeft(), IsApproxEqual(Vector<2>{{-M_SQRT1_2, M_SQRT1_2}}, 1e-12));
        REQUIRE_THAT(bounding.getTopRight(), IsApproxEqual(Vector<2>{{-0.5, std::sqrt(3)/2}}, 1e-12));
    }

    SECTION("only second quarter") {
        RectangularBounding bounding = RectangularBoundingBuilder::forArch({{0, 1}}, 2 * M_PI / 3, 3 * M_PI / 4);

        REQUIRE_THAT(bounding.getBottomLeft(), IsApproxEqual(Vector<2>{{-std::sqrt(3)/2, -M_SQRT1_2}}, 1e-12));
        REQUIRE_THAT(bounding.getTopRight(), IsApproxEqual(Vector<2>{{-M_SQRT1_2, -0.5}}, 1e-12));
    }

    SECTION("only fourth quarter") {
        RectangularBounding bounding = RectangularBoundingBuilder::forArch({{0, 1}}, 5 * M_PI / 3, 7 * M_PI / 4);

        REQUIRE_THAT(bounding.getBottomLeft(), IsApproxEqual(Vector<2>{{M_SQRT1_2, 0.5}}, 1e-12));
        REQUIRE_THAT(bounding.getTopRight(), IsApproxEqual(Vector<2>{{std::sqrt(3)/2, M_SQRT1_2}}, 1e-12));
    }

    SECTION("first and second quarter") {
        RectangularBounding bounding = RectangularBoundingBuilder::forArch({{0, 1}}, M_PI / 6, 3 * M_PI / 4);

        REQUIRE_THAT(bounding.getBottomLeft(), IsApproxEqual(Vector<2>{{-1, -M_SQRT1_2}}, 1e-12));
        REQUIRE_THAT(bounding.getTopRight(), IsApproxEqual(Vector<2>{{-0.5, std::sqrt(3)/2}}, 1e-12));
    }

    SECTION("third and fourth quarter") {
        RectangularBounding bounding = RectangularBoundingBuilder::forArch({{0, 1}}, 5 * M_PI / 4, 5 * M_PI / 3);

        REQUIRE_THAT(bounding.getBottomLeft(), IsApproxEqual(Vector<2>{{M_SQRT1_2, -M_SQRT1_2}}, 1e-12));
        REQUIRE_THAT(bounding.getTopRight(), IsApproxEqual(Vector<2>{{1, 0.5}}, 1e-12));
    }

    SECTION("second, third, fourth and fifth quarter") {
        RectangularBounding bounding = RectangularBoundingBuilder::forArch({{0, 1}}, 2 * M_PI / 3, 9 * M_PI / 4);

        REQUIRE_THAT(bounding.getBottomLeft(), IsApproxEqual(Vector<2>{{-std::sqrt(3)/2, -1}}, 1e-12));
        REQUIRE_THAT(bounding.getTopRight(), IsApproxEqual(Vector<2>{{1, 1}}, 1e-12));
    }
}