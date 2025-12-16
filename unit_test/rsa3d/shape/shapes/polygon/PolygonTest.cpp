//
// Created by pkua on 03.05.19.
//

#include <catch.hpp>
#include "../../../../../rsa3d/shape/shapes/polygon/Polygon.h"

TEST_CASE("Polygon: initClass basic parsing") {
    SECTION("xy 1x1 square vertices") {
        Polygon::initClass("4 xy 0.5 0.5 -0.5 0.5 -0.5 -0.5 0.5 -0.5 4 0 1 2 3 0");

        auto vertexR = Polygon::getVertexR();
        auto vertexTheta = Polygon::getVertexTheta();
        REQUIRE(vertexR.size() == 4);
        REQUIRE(vertexTheta.size() == 4);
        REQUIRE(vertexR[0] == Approx(M_SQRT1_2));
        REQUIRE(vertexR[1] == Approx(M_SQRT1_2));
        REQUIRE(vertexR[2] == Approx(M_SQRT1_2));
        REQUIRE(vertexR[3] == Approx(M_SQRT1_2));
        REQUIRE(vertexTheta[0] == Approx(M_PI/4));
        REQUIRE(vertexTheta[1] == Approx(3*M_PI/4));
        REQUIRE(vertexTheta[2] == Approx(-3*M_PI/4));
        REQUIRE(vertexTheta[3] == Approx(-M_PI/4));
    }

    SECTION("rt 1x1 square vertices") {
        std::ostringstream attr;
        attr << "4 rt " << M_SQRT1_2 << " " << M_PI/4 << " " << M_SQRT1_2 << " " << 3*M_PI/4 << " " << M_SQRT1_2 << " ";
        attr << -3*M_PI/4 << " " << M_SQRT1_2 << " " << -M_PI/4 << " 4 0 1 2 3 0";

        Polygon::initClass(attr.str());

        auto vertexR = Polygon::getVertexR();
        auto vertexTheta = Polygon::getVertexTheta();
        REQUIRE(vertexR.size() == 4);
        REQUIRE(vertexTheta.size() == 4);
        REQUIRE(vertexR[0] == Approx(M_SQRT1_2));
        REQUIRE(vertexR[1] == Approx(M_SQRT1_2));
        REQUIRE(vertexR[2] == Approx(M_SQRT1_2));
        REQUIRE(vertexR[3] == Approx(M_SQRT1_2));
        REQUIRE(vertexTheta[0] == Approx(M_PI/4));
        REQUIRE(vertexTheta[1] == Approx(3*M_PI/4));
        REQUIRE(vertexTheta[2] == Approx(-3*M_PI/4));
        REQUIRE(vertexTheta[3] == Approx(-M_PI/4));
    }

    SECTION("xy 1x1 square segments") {
        Polygon::initClass("4 xy 0.5 0.5 -0.5 0.5 -0.5 -0.5 0.5 -0.5 4 0 1 2 3 0");

        auto segments = Polygon::getSegments();
        auto helperSegments = Polygon::getHelperSegments();
        REQUIRE(segments.size() == 4);
        REQUIRE(segments[0] == std::make_pair((size_t)0, (size_t)1));
        REQUIRE(segments[1] == std::make_pair((size_t)1, (size_t)2));
        REQUIRE(segments[2] == std::make_pair((size_t)2, (size_t)3));
        REQUIRE(segments[3] == std::make_pair((size_t)3, (size_t)0));
        REQUIRE(helperSegments.empty());
    }
}

TEST_CASE("Polygon: volume normalization") {
    Polygon::initClass("4 xy 1.3 1.3 -1.3 1.3 -1.3 -1.3 1.3 -1.3 4 0 1 2 3 0");

    auto vertexR = Polygon::getVertexR();
    auto vertexTheta = Polygon::getVertexTheta();
    REQUIRE(vertexR.size() == 4);
    REQUIRE(vertexTheta.size() == 4);
    REQUIRE(vertexR[0] == Approx(M_SQRT1_2));
    REQUIRE(vertexR[1] == Approx(M_SQRT1_2));
    REQUIRE(vertexR[2] == Approx(M_SQRT1_2));
    REQUIRE(vertexR[3] == Approx(M_SQRT1_2));
    REQUIRE(vertexTheta[0] == Approx(M_PI / 4));
    REQUIRE(vertexTheta[1] == Approx(3 * M_PI / 4));
    REQUIRE(vertexTheta[2] == Approx(-3 * M_PI / 4));
    REQUIRE(vertexTheta[3] == Approx(-M_PI / 4));
}

TEST_CASE("Polygon: automatic volume only for convex") {
    REQUIRE_THROWS_WITH(Polygon::initClass("4 xy 1.3 1.3 -1.3 -1.3 -1.3 1.3 1.3 -1.3 4 0 1 2 3 0"),
                        Catch::Contains("unsupported"));
}

TEST_CASE("Polygon: insphere and circumsphere") {
    SECTION("1 x 1 square, center in x = 0.1, y = 0.2") {
        Polygon::initClass("4 xy 0.6 0.7 -0.4 0.7 -0.4 -0.3 0.6 -0.3 4 0 1 2 3 0");

        REQUIRE(Shape<2, 1>::getInsphereRadius() == Approx(0.5).margin(1e-7));
        REQUIRE(Shape<2, 1>::getCircumsphereRadius() == Approx(M_SQRT1_2).margin(1e-7));
    }

    SECTION("3 x 4 x 5 right triangle, center in right angle vertex") {
        // Use external volume = 1 not to normalize anything
        Polygon::initClass("3 xy  0 0 4 0 0 3  3 0 1 2  0  1");

        REQUIRE(Shape<2, 1>::getInsphereRadius() == Approx(1).margin(1e-7));
        REQUIRE(Shape<2, 1>::getCircumsphereRadius() == Approx(std::sqrt(1*1 + 3*3)).margin(1e-7));
    }

    SECTION("M-like shape with a base") {
        // Use external volume = 1 not to normalize anything
        Polygon::initClass("5 xy  0 0 2 0 2 5 1 2 0 5  5 0 1 2 3 4  0  1");

        REQUIRE(Shape<2, 1>::getInsphereRadius() == Approx(1).margin(1e-7));
        REQUIRE(Shape<2, 1>::getCircumsphereRadius() == Approx(std::sqrt(1*1 + 4*4)).margin(1e-7));
    }
}

TEST_CASE("Polygon: initClass helperSegments") {
    SECTION("manual helper segments parsing") {
        Polygon::initClass("4 xy 0.5 0.5 -0.5 0.5 -0.5 -0.5 0.5 -0.5 4 0 1 2 3 2 3 1 2 0");

        auto helperSegments = Polygon::getHelperSegments();
        REQUIRE(helperSegments.size() == 2);
        REQUIRE(helperSegments[0] == std::make_pair((size_t)3, (size_t)1));
        REQUIRE(helperSegments[1] == std::make_pair((size_t)2, (size_t)0));
    }

    SECTION("starHelper segments additional vertex") {
        Polygon::initClass("4 xy 0.5 0.5 -0.5 0.5 -0.5 -0.5 0.5 -0.5 4 0 1 2 3 starHelperSegments");

        auto vertexR = Polygon::getVertexR();
        auto vertexTheta = Polygon::getVertexTheta();
        REQUIRE(vertexR.size() == 5);
        REQUIRE(vertexTheta.size() == 5);
        REQUIRE(vertexR[4] == 0);
        REQUIRE(vertexTheta[4] == 0);
    }

    SECTION("starHelper segments") {
        Polygon::initClass("4 xy 0.5 0.5 -0.5 0.5 -0.5 -0.5 0.5 -0.5 4 0 1 2 3 starHelperSegments");

        auto helperSegments = Polygon::getHelperSegments();
        REQUIRE(helperSegments.size() == 4);
        REQUIRE(helperSegments[0] == std::make_pair((size_t)0, (size_t)4));
        REQUIRE(helperSegments[1] == std::make_pair((size_t)1, (size_t)4));
        REQUIRE(helperSegments[2] == std::make_pair((size_t)2, (size_t)4));
        REQUIRE(helperSegments[3] == std::make_pair((size_t)3, (size_t)4));
    }
}

TEST_CASE("Polygon: getVolume") {
    SECTION("getVolume should always give 1") {
        Polygon::initClass("4 xy 0.5 0.5 -0.5 0.5 -0.5 -0.5 0.5 -0.5 4 0 1 2 3 0");
        Polygon polygon1;
        REQUIRE(polygon1.getVolume(2) == Approx(1));

        Polygon::initClass("4 xy 0.6 0.7 -0.4 0.7 -0.4 -0.3 0.6 -0.3 4 0 1 2 3 0");
        Polygon polygon2;
        REQUIRE(polygon2.getVolume(2) == Approx(1));

        Polygon::initClass("4 xy 0.6 0.9 -0.5 0.7 -0.4 -0.2 0.6 -0.3 4 0 1 2 3 0");
        Polygon polygon3;
        REQUIRE(polygon3.getVolume(2) == Approx(1));
    }
}