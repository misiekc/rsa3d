//
// Created by pkua on 15.11.2019.
//

#include <catch.hpp>

#include "../../../../../rsa3d/shape/shapes/polydisk/Kmer.h"

TEST_CASE("Kmer: argument parsing exceptions") {
    SECTION("no or too little arguments should throw") {
        REQUIRE_THROWS_WITH(Kmer::initClass(""), Catch::Contains("arguments format"));
        REQUIRE_THROWS_WITH(Kmer::initClass("1"), Catch::Contains("arguments format"));
    }

    SECTION("number of disks >= 2") {
        REQUIRE_THROWS_WITH(Kmer::initClass("1 1"), Catch::Contains(">= 2"));
        REQUIRE_NOTHROW(Kmer::initClass("2 1"));
    }

    SECTION("distance >= 0") {
        REQUIRE_THROWS_WITH(Kmer::initClass("2 -0.1"), Catch::Contains(">= 0"));
        REQUIRE_NOTHROW(Kmer::initClass("2 0"));
    }
}

TEST_CASE("Kmer: preparing Polydisk attributes") {
    SECTION("3 disks, distance 0.5") {
        std::string attr = Kmer::preparePolydiskAttr(3, 0.5);

        // Stripping area
        REQUIRE(attr.length() > 32);
        std::string expected = "3 xy -0.5 0 1 0 0 1 0.5 0 1";
        REQUIRE(attr.substr(0, expected.length()) == expected);
    }

    SECTION("4 disks, distance 0.5") {
        std::string attr = Kmer::preparePolydiskAttr(4, 0.5);

        // Stripping area
        REQUIRE(attr.length() > 32);
        std::string expected = "4 xy -0.75 0 1 -0.25 0 1 0.25 0 1 0.75 0 1";
        REQUIRE(attr.substr(0, expected.length()) == expected);
    }
}

TEST_CASE("Kmer: calculating area") {
    SECTION("2 disks, distance 2, touching") {
        REQUIRE(Kmer::calculateArea(2, 2) == Approx(6.60569426872755));
    }

    SECTION("4 disks, distance 0.5 - overlapping") {
        REQUIRE(Kmer::calculateArea(4, 0.5) == Approx(6.11806287853746));
    }

    SECTION("2 disks, distance 2.5 - separated") {
        REQUIRE(Kmer::calculateArea(2, 2.5) == Approx(7.04471640258879));
    }

    SECTION("2 disks, distance 2sqrt(2) - separated as far as possible") {
        REQUIRE(Kmer::calculateArea(2, 2*M_SQRT2) == Approx(4 + M_PI));
    }
}