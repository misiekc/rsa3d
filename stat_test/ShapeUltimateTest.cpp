//
// Created by pkua on 07.06.19.
//

#include <map>
#include <cstring>
#include <iomanip>

#include "ShapeUltimateTest.h"
#include "ShapeOverlapTest.h"
#include "ShapeStaticSizesTest.h"
#include "ShapePointInsideTest.h"
#include "ShapeBCTest.h"
#include "utils/TestExitCodes.h"

namespace {
    using MainFunc = decltype(&shape_ultitest::main);

    struct TestInfo {
        const std::string name;
        const MainFunc mainFunc;
        const std::vector<std::string> arguments;
        int result = TEST_FAILURE;

        void print(std::ostream &out) const {
            std::string resultDescription;
            switch(result) {
                case TEST_SUCCESS:
                    resultDescription = "SUCCESS";
                    break;
                case TEST_FAILURE:
                    resultDescription = "FAILURE";
                    break;
                case TEST_ERROR:
                    resultDescription = "INTERNAL ERROR";
                    break;
                case TEST_UNSUPPORTED:
                    resultDescription = "UNSUPPORTED";
                    break;
                default:
                    resultDescription = "UNKNOWN";
                    break;
            }

            out << std::setw(16) << std::left << name << " : " << resultDescription << std::endl;
        }
    };

    std::vector<TestInfo> testInfos = {
            {"shape_ovtest", shape_ovtest::main, {"1.7", "100000"}},
            {"shape_pitest", shape_pitest::main, {"1.7", "100000"}},
            {"shape_sizetest", shape_sizetest::main, {"100000"}},
            {"shape_bctest", shape_bctest::main, {"1.7", "100000"}}
    };

    class ArgumentsBuilder {
    private:
        int argc;
        char **argv;

    public:
        explicit ArgumentsBuilder(const std::vector<std::string> &arguments) {
            argc = arguments.size();
            argv = new char*[argc];

            for (int i = 0; i < argc; i++) {
                auto argument = arguments[i];
                argv[i] = new char[argument.length() + 1];
                std::copy(argument.begin(), argument.end(), argv[i]);
                argv[i][argument.length()] = '\0';
            }
        }

        ~ArgumentsBuilder() {
            for (int i = 0; i < argc; i++)
                delete [] argv[i];
            delete [] argv;
        }

        int getArgc() const { return argc; }
        char **getArgv() const { return argv; }
    };
}

int shape_ultitest::main(int argc, char **argv) {
    if (argc < 4) {
        std::cerr << "Usage: ./rsa_test shape_ovtest [particle] [attibutes]" << std::endl;
        return TEST_ERROR;
    }

    std::string particle = argv[2];
    std::string attributes = argv[3];

    for (auto &testInfo : testInfos) {
        std::vector<std::string> arguments = {argv[0], testInfo.name, particle, attributes};
        std::copy(testInfo.arguments.begin(), testInfo.arguments.end(), std::back_inserter(arguments));
        ArgumentsBuilder builder(arguments);

        std::cout << std::endl;
        std::cout << "************************************************************" << std::endl;
        std::cout << "* Starting " << testInfo.name << "..." << std::endl;
        std::cout << "************************************************************" << std::endl << std::endl;

        testInfo.result = testInfo.mainFunc(builder.getArgc(), builder.getArgv());
    }

    std::cout << std::endl;
    std::cout << "************************************************************" << std::endl;
    std::cout << "* Results summary" << std::endl;
    std::cout << "************************************************************" << std::endl << std::endl;
    for (const auto& testInfo : testInfos)
        testInfo.print(std::cout);

    std::size_t numFails = std::count_if(testInfos.begin(), testInfos.end(), [](auto info) {
        return info.result == TEST_FAILURE;
    });

    std::size_t numError = std::count_if(testInfos.begin(), testInfos.end(), [](auto info) {
        return info.result == TEST_ERROR;
    });

    if (numError != 0)
        return TEST_ERROR;
    else if (numFails != 0)
        return TEST_FAILURE;
    else
        return  TEST_SUCCESS;
}

