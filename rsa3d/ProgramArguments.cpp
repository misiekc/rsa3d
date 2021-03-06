//
// Created by pkua on 07.04.19.
//

#include <fstream>
#include <iostream>

#include "ProgramArguments.h"
#include "utils/Assertions.h"

namespace {
    /* The obvious implementation of reading data from file. It is extracted only for testability (mocking of reading
     * from file. */
    class ParametersFileReaderImpl : public ParametersFileReader {
    public:
        ~ParametersFileReaderImpl() override = default;

        /* It just opens file with name filename and stores it into output stream. */
        void readFromFile(const std::string &filename, std::ostream &output) override {
            std::ifstream inputFile(filename);
            if (!inputFile)
                throw InvalidArgumentsException("Cannot open input file " + filename + " to read. Aborting.");

            output << inputFile.rdbuf() << std::endl;
        }
    };
}

ProgramArguments::ProgramArguments(int argc, char **argv) {
    ParametersFileReaderImpl pfr;
    *this = ProgramArguments(argc, argv, pfr, std::cin);
}

ProgramArguments::ProgramArguments(int argc, char **argv, ParametersFileReader &pfr, std::istream &inputStream) {
    Expects(argc >= 1);     // valid program arguments have at least 1 arg - cmd; assertion error otherwise
    this->cmd = argv[0];

    std::stringstream input;
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-f") {
            if (++i == argc)
                throw InvalidArgumentsException("Expected input file after -f. Aborting.");
            pfr.readFromFile(argv[i], input);
            input << std::endl;
        } else if (arg == "-i") {
            input << inputStream.rdbuf() << std::endl;
        } else if (startsWithMinus(arg)) {
            input << arg.substr(1) << std::endl;
        } else {
            this->positionalArguments.push_back(arg);
        }
    }

    this->fetchModeArgument();
    this->parameters = Parameters(input);
}

bool ProgramArguments::startsWithMinus(const std::string &arg) const {
    return arg.length() > 1 && arg[0] == '-';
}

void ProgramArguments::fetchModeArgument() {
    // First positional argument, after stripping all extra stuff, should be the mode parameter
    if (this->positionalArguments.empty())
        throw InvalidArgumentsException(this->formatUsage(""));

    this->mode = this->positionalArguments[0];
    this->positionalArguments.erase(this->positionalArguments.begin());

    // If the mode is help, then we want to flag, that help is wanted and mode is actually the argument after help
    // (or empty string, if none). Then the class is in the state "help wanted for this mode"
    if (this->mode == "help") {
        this->helpRequested = true;
        if (this->positionalArguments.empty()) {
            this->mode = "";
        } else {
            this->mode = this->positionalArguments[0];
            this->positionalArguments.erase(this->positionalArguments.begin());
        }
    }
}

std::string ProgramArguments::formatUsage(const std::string &additionalArgs) const {
    std::ostringstream out;
    out << "Usage: " << this->cmd << " ";

    // If mode is unknown, we use a placeholder and put "mode specific arguments" instead of additionalArgs
    if (this->mode.empty())
        out << "[mode] (mode specific arguments)";
    else
        out << this->mode << " " << additionalArgs;
    out << " (flag parameters anywhere)";

    return out.str();
}
