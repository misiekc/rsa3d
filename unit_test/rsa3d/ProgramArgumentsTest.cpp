//
// Created by pkua on 07.04.19.
//

#include <catch.hpp>
#include "../../rsa3d/ProgramArguments.h"

TEST_CASE("ProgramArguments: simple") {
    SECTION("No cmd and param mode should throw") {
        char **argv1 = nullptr;
        char *argv2[] = {(char*)"cmd"};

        REQUIRE_THROWS(ProgramArguments{0, argv1});
        REQUIRE_THROWS(ProgramArguments{1, argv2});
    }

    SECTION("Only mode") {
        char *argv[] = {(char*)"cmd", (char*)"mode"};

        ProgramArguments pa{2, argv};

        REQUIRE(pa.getParameters() == Parameters{});    // default
        REQUIRE(pa.getCmd() == "cmd");
        REQUIRE(pa.getMode() == "mode");
    }
}

TEST_CASE("ProgramArguments: command line parameters") {
    char *argv[] = {(char*)"cmd", (char*)"mode", (char*)"-particleType=Ellipse"};
    Parameters expectedParameters;
    expectedParameters.particleType = "Ellipse";

    ProgramArguments pa{3, argv};

    REQUIRE(pa.getParameters() == expectedParameters);
    REQUIRE(pa.getCmd() == "cmd");
    REQUIRE(pa.getMode() == "mode");
}

TEST_CASE("ProgramArguments: positional arguments") {
    char *argv[] = {(char*)"cmd", (char*)"mode", (char*)"pos1", (char*)"-particleType=Ellipse", (char*)"pos2"};
    Parameters expectedParameters;
    expectedParameters.particleType = "Ellipse";

    ProgramArguments pa{5, argv};

    auto posArgs = pa.getPositionalArguments();
    REQUIRE(pa.getParameters() == expectedParameters);
    REQUIRE(pa.getCmd() == "cmd");
    REQUIRE(pa.getMode() == "mode");
    REQUIRE(posArgs.size() == 2);
    REQUIRE(posArgs[0] == "pos1");
    REQUIRE(posArgs[1] == "pos2");
}