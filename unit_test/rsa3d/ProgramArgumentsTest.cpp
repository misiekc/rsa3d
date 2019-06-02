//
// Created by pkua on 07.04.19.
//

#include <catch.hpp>
#include "../../rsa3d/ProgramArguments.h"

using namespace Catch::Matchers;

TEST_CASE("ProgramArguments: simple") {
    SECTION("No cmd and param mode should throw") {
        char **argv1 = nullptr;
        char *argv2[] = {(char*)"cmd"};
        char *argv3[] = {(char*)"cmd", (char*)"-particleType=Rectangle"};

        REQUIRE_THROWS(ProgramArguments{0, argv1});
        REQUIRE_THROWS(ProgramArguments{1, argv2});
        REQUIRE_THROWS(ProgramArguments{2, argv3});
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
    SECTION("no overriding") {
        char *argv[] = {(char *) "cmd", (char *) "mode", (char *) "-particleType=Ellipse"};
        Parameters expectedParameters;
        expectedParameters.particleType = "Ellipse";

        ProgramArguments pa{3, argv};

        REQUIRE(pa.getParameters() == expectedParameters);
        REQUIRE(pa.getCmd() == "cmd");
        REQUIRE(pa.getMode() == "mode");
    }

    SECTION("overriding") {
        char *argv[] = {(char *) "cmd", (char *) "mode", (char *) "-particleType=Ellipse",
                        (char *) "-particleType=Polygon"};
        Parameters expectedParameters;
        expectedParameters.particleType = "Polygon";

        ProgramArguments pa{4, argv};

        REQUIRE(pa.getParameters() == expectedParameters);
        REQUIRE(pa.getCmd() == "cmd");
        REQUIRE(pa.getMode() == "mode");
    }
}

class MockParametersFileReader : public ParametersFileReader {
public:
    void readFromFile(const std::string &filename, std::ostream &output) override {
        REQUIRE(filename == "file");
        output << "particleType=Ellipse" << std::endl << "particleAttributes=1.5";
    }
};

class MockParametersFileReaderNoOverriding : public ParametersFileReader {
public:
    void readFromFile(const std::string &filename, std::ostream &output) override {
        REQUIRE((filename == "file" || filename == "no_ov"));
        if (filename == "file")
            output << "particleType=Ellipse" << std::endl << "particleAttributes=1.5";
        else if (filename == "no_ov")
            output << "ompThreads=10";
    }
};

class MockParametersFileReaderOverriding : public ParametersFileReader {
public:
    void readFromFile(const std::string &filename, std::ostream &output) override {
        REQUIRE((filename == "file" || filename == "ov"));
        if (filename == "file")
            output << "particleType=Ellipse" << std::endl << "particleAttributes=1.5";
        else if (filename == "ov")
            output << "particleType=Polygon";
    }
};

TEST_CASE("ProgramArguments: file parameters") {
    SECTION("valid") {
        char *argv[] = {(char *) "cmd", (char *) "mode", (char *) "-f", (char *) "file"};
        MockParametersFileReader pfr;
        std::istringstream stream;
        Parameters expectedParameters;
        expectedParameters.particleType = "Ellipse";
        expectedParameters.particleAttributes = "1.5";

        ProgramArguments pa{4, argv, pfr, stream};

        REQUIRE(pa.getParameters() == expectedParameters);
        REQUIRE(pa.getCmd() == "cmd");
        REQUIRE(pa.getMode() == "mode");
    }

    SECTION("no file should throw") {
        char *argv[] = {(char *) "cmd", (char *) "mode", (char *) "-f"};

        REQUIRE_THROWS_WITH(ProgramArguments(3, argv), Contains("Expected input file after -f"));
    }

    SECTION("two files not overriding") {
        char *argv[] = {(char *) "cmd", (char *) "mode", (char *) "-f", (char *) "file", (char *) "-f",
                        (char *) "no_ov"};
        MockParametersFileReaderNoOverriding pfr;
        std::istringstream stream;
        Parameters expectedParameters;
        expectedParameters.particleType = "Ellipse";
        expectedParameters.particleAttributes = "1.5";
        expectedParameters.ompThreads = 10;

        ProgramArguments pa{6, argv, pfr, stream};

        REQUIRE(pa.getParameters() == expectedParameters);
        REQUIRE(pa.getCmd() == "cmd");
        REQUIRE(pa.getMode() == "mode");
    }

    SECTION("two files overriding") {
        char *argv[] = {(char *) "cmd", (char *) "mode", (char *) "-f", (char *) "file", (char *) "-f",
                        (char *) "ov"};
        MockParametersFileReaderOverriding pfr;
        std::istringstream stream;
        Parameters expectedParameters;
        expectedParameters.particleType = "Polygon";
        expectedParameters.particleAttributes = "1.5";

        ProgramArguments pa{6, argv, pfr, stream};

        REQUIRE(pa.getParameters() == expectedParameters);
        REQUIRE(pa.getCmd() == "cmd");
        REQUIRE(pa.getMode() == "mode");
    }
}

TEST_CASE("ProgramArguments: stdin parameters") {
    char *argv[] = {(char *) "cmd", (char *) "mode", (char *) "-i"};
    MockParametersFileReader pfr;
    std::istringstream stream("particleType=Ellipse");
    Parameters expectedParameters;
    expectedParameters.particleType = "Ellipse";

    ProgramArguments pa{3, argv, pfr, stream};

    REQUIRE(pa.getParameters() == expectedParameters);
    REQUIRE(pa.getCmd() == "cmd");
    REQUIRE(pa.getMode() == "mode");
}

TEST_CASE("ProgramArguments: mixed parameters") {
    SECTION("no overriding") {
        char *argv[] = {(char *) "cmd", (char *) "mode", (char *) "-f", (char *) "file", (char*) "-ompThreads=10",
                        (char *) "-i"};
        MockParametersFileReader pfr;
        std::istringstream stream("boundaryConditions=free");
        Parameters expectedParameters;
        expectedParameters.particleType = "Ellipse";
        expectedParameters.particleAttributes = "1.5";
        expectedParameters.ompThreads = 10;
        expectedParameters.boundaryConditions = "free";

        ProgramArguments pa{6, argv, pfr, stream};

        REQUIRE(pa.getParameters() == expectedParameters);
        REQUIRE(pa.getCmd() == "cmd");
        REQUIRE(pa.getMode() == "mode");
    }

    SECTION("no overriding + mode parameter at the end") {
        char *argv[] = {(char *) "cmd", (char *) "-f", (char *) "file", (char*) "-ompThreads=10", (char *) "-i",
                        (char *) "mode"};
        MockParametersFileReader pfr;
        std::istringstream stream("boundaryConditions=free");
        Parameters expectedParameters;
        expectedParameters.particleType = "Ellipse";
        expectedParameters.particleAttributes = "1.5";
        expectedParameters.ompThreads = 10;
        expectedParameters.boundaryConditions = "free";

        ProgramArguments pa{6, argv, pfr, stream};

        REQUIRE(pa.getParameters() == expectedParameters);
        REQUIRE(pa.getCmd() == "cmd");
        REQUIRE(pa.getMode() == "mode");
    }

    SECTION("cmd overriding input+file") {
        char *argv[] = {(char *) "cmd", (char *) "mode", (char*) "-i", (char *) "-f", (char *) "file",
                        (char *) "-ompThreads=6", (char *) "-particleAttributes=2.0"};
        MockParametersFileReader pfr;
        std::istringstream stream("ompThreads=5");
        Parameters expectedParameters;
        expectedParameters.particleType = "Ellipse";
        expectedParameters.particleAttributes = "2.0";
        expectedParameters.ompThreads = 6;

        ProgramArguments pa{7, argv, pfr, stream};

        REQUIRE(pa.getParameters() == expectedParameters);
        REQUIRE(pa.getCmd() == "cmd");
        REQUIRE(pa.getMode() == "mode");
    }

    SECTION("file overriding cmd+input") {
        char *argv[] = {(char *) "cmd", (char *) "mode", (char *) "-particleType=Polygon", (char *) "-i", (char *) "-f",
                        (char *) "file"};
        MockParametersFileReader pfr;
        std::istringstream stream("particleAttributes=4.5");
        Parameters expectedParameters;
        expectedParameters.particleType = "Ellipse";
        expectedParameters.particleAttributes = "1.5";

        ProgramArguments pa{6, argv, pfr, stream};

        REQUIRE(pa.getParameters() == expectedParameters);
        REQUIRE(pa.getCmd() == "cmd");
        REQUIRE(pa.getMode() == "mode");
    }

    SECTION("input overriding file+cmd") {
        char *argv[] = {(char *) "cmd", (char *) "mode", (char *) "-f", (char *) "file",
                        (char *) "-ompThreads=10", (char *) "-i"};
        MockParametersFileReader pfr;
        std::stringstream stream;
        stream << "particleType=Spherocylinder2D" << std::endl << "ompThreads=9";
        Parameters expectedParameters;
        expectedParameters.particleType = "Spherocylinder2D";
        expectedParameters.particleAttributes = "1.5";
        expectedParameters.ompThreads = 9;

        ProgramArguments pa{6, argv, pfr, stream};

        REQUIRE(pa.getParameters() == expectedParameters);
        REQUIRE(pa.getCmd() == "cmd");
        REQUIRE(pa.getMode() == "mode");
    }
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