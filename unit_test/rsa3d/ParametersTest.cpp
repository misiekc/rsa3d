//
// Created by PKua on 30.09.18.
//

//#####################//
// ENABLED FOR rsa.2.0 //
//#####################//

#if RSA_SPATIAL_DIMENSION == 2 && RSA_ANGULAR_DIMENSION == 0

#include <catch.hpp>
#include "../../rsa3d/Parameters.h"


TEST_CASE("Parameters: valid") {
    SECTION("default construction should not throws") {
        REQUIRE_NOTHROW(Parameters{});
    }

    SECTION("all parameters != default") {
        std::istringstream str{"maxTriesWithoutSuccess = 1\n"
                               "maxVoxels = 2\n"
                               "requestedAngularVoxelSize = 3.5\n"
                               "minDx = 4.5\n"
                               "maxTime = 5.5\n"
                               "split = 6\n"
                               "surfaceVolume = 2.25\n"
                               "storePackings = false\n"
                               "modifiedRSA = true\n"
                               "thresholdDistance = 7.5\n"
                               "boundaryConditions = free\n"
                               "particleType = Particle\n"
                               "particleAttributes = attr\n"
                               "from = 8\n"
                               "collectors = 9\n"
                               "generatorProcesses = 10\n"
                               "ompThreads = 11\n"};

        Parameters parameters{str};

        REQUIRE(parameters.maxTriesWithoutSuccess == 1);
        REQUIRE(parameters.maxVoxels == 2);
        REQUIRE(parameters.requestedAngularVoxelSize == Approx(3.5));
        REQUIRE(parameters.minDx == Approx(4.5));
        REQUIRE(parameters.maxTime == Approx(5.5));
        REQUIRE(parameters.split == 6);
        REQUIRE(parameters.surfaceSize == Approx(1.5));
        REQUIRE(parameters.storePackings == false);
        REQUIRE(parameters.modifiedRSA == true);
        REQUIRE(parameters.thresholdDistance == Approx(7.5));
        REQUIRE(parameters.boundaryConditions == "free");
        REQUIRE(parameters.particleType == "Particle");
        REQUIRE(parameters.particleAttributes == "attr");
        REQUIRE(parameters.from == 8);
        REQUIRE(parameters.collectors == 9);
        REQUIRE(parameters.generatorProcesses == 10);
        REQUIRE(parameters.ompThreads == 11);
    }
    
    SECTION("some parameters, should not throw") {
        std::istringstream str{"maxTriesWithoutSuccess = 1\n"
                               "minDx = 4.5\n"
                               "maxTime = 5.5\n"
                               "split = 6\n"
                               "generatorProcesses = 10\n"
                               "ompThreads = 11\n"};

        REQUIRE_NOTHROW(Parameters{str});
    }
    
    SECTION("some 0s allowed, should not throw") {
        std::istringstream str{"minDx = 0\n"
                               "particleAttributes =\n"
                               "from = 0\n"
                               "thresholdDistance = 0\n"};

        REQUIRE_NOTHROW(Parameters{str});   
    }
}

TEST_CASE("Parameters: invalid") {
    SECTION("maxTriesWithoutSuccess > 0") {
        std::istringstream str1{"maxTriesWithoutSuccess = 0"};
        std::istringstream str2{"maxTriesWithoutSuccess = -1"};

        REQUIRE_THROWS(Parameters{str1});
        REQUIRE_THROWS(Parameters{str2});
    };

    SECTION("maxVoxels >= 0") {
        std::istringstream str2{"maxVoxels = -1"};

        REQUIRE_THROWS(Parameters{str2});
    }

    SECTION("requestedAngularVoxelSize > 0") {
        std::istringstream str1{"requestedAngularVoxelSize = 0"};
        std::istringstream str2{"requestedAngularVoxelSize = -0.5"};

        REQUIRE_THROWS(Parameters{str1});
        REQUIRE_THROWS(Parameters{str2});
    }

    SECTION("minDx >= 0.0") {
        std::istringstream str{"minDx = -0.5"};

        REQUIRE_THROWS(Parameters{str});
    }

    SECTION("from >= 0") {
        std::istringstream str{"from = -1"};

        REQUIRE_THROWS(Parameters{str});
    }

    SECTION("collectors > 0") {
        std::istringstream str1{"collectors = 0"};
        std::istringstream str2{"collectors = -1"};

        REQUIRE_THROWS(Parameters{str1});
        REQUIRE_THROWS(Parameters{str2});
    }

    SECTION("maxTime > 0") {
        std::istringstream str1{"maxTime = 0"};
        std::istringstream str2{"maxTime = -1"};

        REQUIRE_THROWS(Parameters{str1});
        REQUIRE_THROWS(Parameters{str2});
    }

    SECTION("split > 0") {
        std::istringstream str1{"split = 0"};
        std::istringstream str2{"split = -1"};

        REQUIRE_THROWS(Parameters{str1});
        REQUIRE_THROWS(Parameters{str2});
    }

    SECTION("!std::isnan(surfaceSize) && surfaceSize > 0") {
        std::istringstream str1{"surfaceVolume = 0"};
        std::istringstream str2{"surfaceVolume = -1"};

        REQUIRE_THROWS(Parameters{str1});
        REQUIRE_THROWS(Parameters{str2});
    }

    SECTION("thresholdDistance >= 0.0") {
        std::istringstream str{"thresholdDistance = -0.5"};

        REQUIRE_THROWS(Parameters{str});
    }

    SECTION("generatorProcesses > 0") {
        std::istringstream str1{"generatorProcesses = 0"};
        std::istringstream str2{"generatorProcesses = -1"};

        REQUIRE_THROWS(Parameters{str1});
        REQUIRE_THROWS(Parameters{str2});
    }

    SECTION("ompThreads > 0") {
        std::istringstream str1{"ompThreads = 0"};
        std::istringstream str2{"ompThreads = -1"};

        REQUIRE_THROWS(Parameters{str1});
        REQUIRE_THROWS(Parameters{str2});
    }

    SECTION("!particleType.empty()") {
        std::istringstream str{"particleType ="};

        REQUIRE_THROWS(Parameters{str});
    }

}

#endif