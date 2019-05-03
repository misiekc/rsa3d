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
        REQUIRE(segments[0] == std::make_pair(0ul, 1ul));
        REQUIRE(segments[1] == std::make_pair(1ul, 2ul));
        REQUIRE(segments[2] == std::make_pair(2ul, 3ul));
        REQUIRE(segments[3] == std::make_pair(3ul, 0ul));
        REQUIRE(helperSegments.empty());
    }
}

TEST_CASE("Polygon: initClass calculations") {
    SECTION("volume normalization") {
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

    SECTION("inscribed and circumscribed sphere radius") {
        // 1 x 1 square, center in x = 0.1, y = 0.2
        Polygon::initClass("4 xy 0.6 0.7 -0.4 0.7 -0.4 -0.3 0.6 -0.3 4 0 1 2 3 0");

        REQUIRE(Polygon::getInscribedCircleRadius() == Approx(0.3));
        REQUIRE(Polygon::getCircumscribedCircleRadius() == Approx(std::sqrt(0.6*0.6 + 0.7*0.7)));
    }
}

TEST_CASE("Polygon: initClass helperSegments") {
    SECTION("manual helper segments parsing") {
        Polygon::initClass("4 xy 0.5 0.5 -0.5 0.5 -0.5 -0.5 0.5 -0.5 4 0 1 2 3 2 3 1 2 0");

        auto helperSegments = Polygon::getHelperSegments();
        REQUIRE(helperSegments.size() == 2);
        REQUIRE(helperSegments[0] == std::make_pair(3ul, 1ul));
        REQUIRE(helperSegments[1] == std::make_pair(2ul, 0ul));
    }

    SECTION("manual helper segments shouldn't center polygon") {
        Polygon::initClass("4 xy 0 0 1 0 1 1 0 1 4 0 1 2 3 2 3 1 2 0");

        auto vertexR = Polygon::getVertexR();
        REQUIRE(vertexR[0] == Approx(0));
        REQUIRE(vertexR[1] == Approx(1));
        REQUIRE(vertexR[2] == Approx(M_SQRT2));
        REQUIRE(vertexR[3] == Approx(1));
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
        REQUIRE(helperSegments[0] == std::make_pair(0ul, 4ul));
        REQUIRE(helperSegments[1] == std::make_pair(1ul, 4ul));
        REQUIRE(helperSegments[2] == std::make_pair(2ul, 4ul));
        REQUIRE(helperSegments[3] == std::make_pair(3ul, 4ul));
    }

    SECTION("star helper segments should center polygon") {
        Polygon::initClass("4 xy 0 0 1 0 1 1 0 1 4 0 1 2 3 starHelperSegments");

        auto vertexR = Polygon::getVertexR();
        REQUIRE(vertexR[0] == Approx(M_SQRT1_2));
        REQUIRE(vertexR[1] == Approx(M_SQRT1_2));
        REQUIRE(vertexR[2] == Approx(M_SQRT1_2));
        REQUIRE(vertexR[3] == Approx(M_SQRT1_2));
    }
}

TEST_CASE("Polygon: getVolume") {
    SECTION("getVolume should always give 1") {
        Polygon::initClass("4 xy 0.5 0.5 -0.5 0.5 -0.5 -0.5 0.5 -0.5 4 0 1 2 3 0");
        Polygon polygon1;
        REQUIRE(polygon1.getVolume() == Approx(1));

        Polygon::initClass("4 xy 0.6 0.7 -0.4 0.7 -0.4 -0.3 0.6 -0.3 4 0 1 2 3 0");
        Polygon polygon2;
        REQUIRE(polygon2.getVolume() == Approx(1));

        Polygon::initClass("4 xy 0.6 0.9 -0.5 0.7 -0.4 -0.2 0.6 -0.3 4 0 1 2 3 0");
        Polygon polygon3;
        REQUIRE(polygon3.getVolume() == Approx(1));
    }
}