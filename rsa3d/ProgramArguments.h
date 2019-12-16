//
// Created by pkua on 07.04.19.
//

#ifndef RSA3D_PROGRAMARGUMENTS_H
#define RSA3D_PROGRAMARGUMENTS_H

#include <string>
#include <utility>
#include <vector>

#include "Parameters.h"


/* Class needed for mocking reading parameters from file in tests */
class ParametersFileReader {
public:
    virtual ~ParametersFileReader() = default;

    virtual void readFromFile(const std::string &filename, std::ostream &output) = 0;
};

class InvalidArgumentsException : public std::runtime_error {
public:
    explicit InvalidArgumentsException(const std::string &msg) : runtime_error{msg} {}
};

class ArgumentsHelpRequest : public std::exception {
private:
    std::string mode;

public:
    explicit ArgumentsHelpRequest(std::string mode) : mode{std::move(mode)} { }

    const std::string &getMode() const { return mode; }
};

class ProgramArguments {
private:
    Parameters parameters;
    std::string cmd;
    std::string mode;
    std::vector<std::string> positionalArguments;

    void fetchModeArgument();
    bool startsWithMinus(const std::string &arg) const;

public:
    ProgramArguments(int argc, char **argv, ParametersFileReader &pfr, std::istream &inputStream);
    ProgramArguments(int argc, char **argv);

    const std::string &getCmd() const { return cmd; }
    const Parameters &getParameters() const { return parameters; }
    const std::string &getMode() const { return mode; }
    const std::vector<std::string> &getPositionalArguments() const { return positionalArguments; }

    std::string formatUsage(const std::string &additionalArgs) const;
};


#endif //RSA3D_PROGRAMARGUMENTS_H
