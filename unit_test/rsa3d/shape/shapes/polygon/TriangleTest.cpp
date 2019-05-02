#include <catch.hpp>
#include "../../../../../rsa3d/shape/shapes/polygon/Triangle.h"

using namespace Catch;

TEST_CASE("Triangle: calculation") {

    SECTION("5 x 13 x 16 triangle output") {
        std::istringstream attr(Triangle::preparePolygonAttributes("25 29 36"));

        // expected output: 3 xy (-18) 0 (18) 0 (-3) (20) 3 0 1 2 0 starHelperSegments  ;  things in ( ) may be approx
        int numV;
        std::string coordType, theRest;
        double v[6];
        attr >> numV >> coordType >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5];
        std::getline(attr, theRest);
        REQUIRE(numV == 3);
        REQUIRE(coordType == "xy");
        REQUIRE(v[0] == Approx(-18));
        REQUIRE(v[1] == 0);
        REQUIRE(v[2] == Approx(18));
        REQUIRE(v[3] == 0);
        REQUIRE(v[4] == Approx(3));
        REQUIRE(v[5] == Approx(20));
        REQUIRE(theRest == " 3 0 1 2 0 starHelperSegments");
    }
}

TEST_CASE("Triangle: normalization") {

    SECTION("a <= b <= c order") {
        std::string attr1 = Triangle::preparePolygonAttributes("3 4 5");
        std::string attr2 = Triangle::preparePolygonAttributes("4 5 3");
        std::string attr3 = Triangle::preparePolygonAttributes("5 3 4");
        std::string attr4 = Triangle::preparePolygonAttributes("5 4 3");
        std::string attr5 = Triangle::preparePolygonAttributes("4 3 5");
        std::string attr6 = Triangle::preparePolygonAttributes("3 5 4");

        REQUIRE(attr1 == attr2);
        REQUIRE(attr2 == attr3);
        REQUIRE(attr3 == attr4);
        REQUIRE(attr4 == attr5);
        REQUIRE(attr5 == attr6);
    }
}

TEST_CASE("Triangle: validation") {
    SECTION("malformed") {
        REQUIRE_THROWS_WITH(Triangle::preparePolygonAttributes("1 2"), Contains("3 doubles"));
        REQUIRE_THROWS_WITH(Triangle::preparePolygonAttributes("a b c"), Contains("3 doubles"));
    }

    SECTION("sides must be > 0") {
        REQUIRE_THROWS_WITH(Triangle::preparePolygonAttributes("-1 1 2"), Contains("positive"));
        REQUIRE_THROWS_WITH(Triangle::preparePolygonAttributes("0 1 2"), Contains("positive"));
        REQUIRE_THROWS_WITH(Triangle::preparePolygonAttributes("1 -1 2"), Contains("positive"));
        REQUIRE_THROWS_WITH(Triangle::preparePolygonAttributes("1 0 2"), Contains("positive"));
        REQUIRE_THROWS_WITH(Triangle::preparePolygonAttributes("1 1 -1"), Contains("positive"));
        REQUIRE_THROWS_WITH(Triangle::preparePolygonAttributes("1 1 0"), Contains("positive"));
    }

    SECTION("sides must not violate triangle equality") {
        REQUIRE_THROWS_WITH(Triangle::preparePolygonAttributes("1 2 3"), Contains("triangle inequality", CaseSensitive::No));
    }
}