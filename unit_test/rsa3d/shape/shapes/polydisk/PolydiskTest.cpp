//
// Created by pkua on 16.09.2019.
//

#include <catch.hpp>
#include "../../../../../rsa3d/shape/shapes/polydisk/Polydisk.h"

class ApproxAngle {
private:
    double angle{};

public:
    ApproxAngle() = default;
    ApproxAngle(double angle) : angle{std::fmod(std::fmod(angle, 2*M_PI) + 2*M_PI, 2*M_PI)} { }

    friend bool operator==(ApproxAngle lhs, ApproxAngle rhs);
};

bool operator==(ApproxAngle lhs, ApproxAngle rhs) {
    return lhs.angle == Approx(rhs.angle);
}

TEST_CASE("Polydisk: basic arguments parsing") {
    SECTION("xy: big circle in centre, 4 small left, right, top, bottom") {
        Polydisk::initClass("5 xy  0 0 2  3 0 1  0 3 1  -3 0 1  0 -3 1  1");

        auto diskCentreR = Polydisk::getDiskCentreR();
        auto diskCentreTheta = Polydisk::getDiskCentreTheta();
        auto diskR = Polydisk::getDiskR();
        REQUIRE(diskCentreR[0] == Approx(0));
        REQUIRE(diskR[0] == Approx(2));

        REQUIRE(diskCentreR[1] == Approx(3));
        REQUIRE(diskCentreTheta[1] == ApproxAngle(0));
        REQUIRE(diskR[1] == Approx(1));

        REQUIRE(diskCentreR[2] == Approx(3));
        REQUIRE(diskCentreTheta[2] == ApproxAngle(M_PI/2));
        REQUIRE(diskR[2] == Approx(1));

        REQUIRE(diskCentreR[3] == Approx(3));
        REQUIRE(diskCentreTheta[3] == ApproxAngle(M_PI));
        REQUIRE(diskR[3] == Approx(1));

        REQUIRE(diskCentreR[4] == Approx(3));
        REQUIRE(diskCentreTheta[4] == ApproxAngle(3*M_PI/2));
        REQUIRE(diskR[4] == Approx(1));
    }

    SECTION("rt: big circle in centre, 4 small approximately left, right, top, bottom") {
        Polydisk::initClass("5 rt  0 0 2  3 0 1  3 1.57 1  3 3.14 1  3 4.71 1  1");

        auto diskCentreR = Polydisk::getDiskCentreR();
        auto diskCentreTheta = Polydisk::getDiskCentreTheta();
        auto diskR = Polydisk::getDiskR();
        REQUIRE(diskCentreR[0] == Approx(0));
        REQUIRE(diskR[0] == Approx(2));

        REQUIRE(diskCentreR[1] == Approx(3));
        REQUIRE(diskCentreTheta[1] == ApproxAngle(0));
        REQUIRE(diskR[1] == Approx(1));

        REQUIRE(diskCentreR[2] == Approx(3));
        REQUIRE(diskCentreTheta[2] == ApproxAngle(1.57));
        REQUIRE(diskR[2] == Approx(1));

        REQUIRE(diskCentreR[3] == Approx(3));
        REQUIRE(diskCentreTheta[3] == ApproxAngle(3.14));
        REQUIRE(diskR[3] == Approx(1));

        REQUIRE(diskCentreR[4] == Approx(3));
        REQUIRE(diskCentreTheta[4] == ApproxAngle(4.71));
        REQUIRE(diskR[4] == Approx(1));
    }
}

TEST_CASE("Polydisk: area scaling") {
    SECTION("simple: small and big circle, scale 2 times (area 4 times)") {
        Polydisk::initClass("2 xy  0 0 4  6 0 2  4");

        auto diskCentreR = Polydisk::getDiskCentreR();
        auto diskCentreTheta = Polydisk::getDiskCentreTheta();
        auto diskR = Polydisk::getDiskR();

        REQUIRE(diskCentreR[0] == Approx(0));
        REQUIRE(diskR[0] == Approx(2));

        REQUIRE(diskCentreR[1] == Approx(3));
        REQUIRE(diskCentreTheta[1] == ApproxAngle(0));
        REQUIRE(diskR[1] == Approx(1));
    }
}

TEST_CASE("Polydisk: centering based on largest circle") {
    SECTION("simple: small and big circle, 0 degrees") {
        Polydisk::initClass("2 xy  -1 0 2  2 0 1  1");

        auto diskCentreR = Polydisk::getDiskCentreR();
        auto diskCentreTheta = Polydisk::getDiskCentreTheta();
        auto diskR = Polydisk::getDiskR();

        REQUIRE(diskCentreR[0] == Approx(0));
        REQUIRE(diskR[0] == Approx(2));

        REQUIRE(diskCentreR[1] == Approx(3));
        REQUIRE(diskCentreTheta[1] == ApproxAngle(0));
        REQUIRE(diskR[1] == Approx(1));
    }

    SECTION("simple: small and big circle, 45 degrees") {
        Polydisk::initClass("2 xy  -1 -1 2  2 2 1  1");

        auto diskCentreR = Polydisk::getDiskCentreR();
        auto diskCentreTheta = Polydisk::getDiskCentreTheta();
        auto diskR = Polydisk::getDiskR();

        REQUIRE(diskCentreR[0] == Approx(0));
        REQUIRE(diskR[0] == Approx(2));

        REQUIRE(diskCentreR[1] == Approx(3*M_SQRT2));
        REQUIRE(diskCentreTheta[1] == ApproxAngle(M_PI/4));
        REQUIRE(diskR[1] == Approx(1));
    }

    SECTION("small and 2 big circle, closer one to choose") {
        Polydisk::initClass("3 xy  -2 0 2  2 0 2  5 0 1  1");

        auto diskCentreR = Polydisk::getDiskCentreR();
        auto diskCentreTheta = Polydisk::getDiskCentreTheta();
        auto diskR = Polydisk::getDiskR();

        REQUIRE(diskCentreR[0] == Approx(4));
        REQUIRE(diskCentreTheta[0] == ApproxAngle(M_PI));
        REQUIRE(diskR[0] == Approx(2));

        REQUIRE(diskCentreR[1] == Approx(0));
        REQUIRE(diskR[1] == Approx(2));

        REQUIRE(diskCentreR[2] == Approx(3));
        REQUIRE(diskCentreTheta[2] == ApproxAngle(0));
        REQUIRE(diskR[2] == Approx(1));
    }
}