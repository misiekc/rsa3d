//
// Created by pkua on 15.05.19.
//

#ifndef RSA3D_ASSERTIONS_H
#define RSA3D_ASSERTIONS_H

#include <stdexcept>

// Cpp Core Guidelines-style assertions for design by contract
// https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#i6-prefer-expects-for-expressing-preconditions

// Preconditions check (argument validation)
#define Expects(cond) if (!(cond)) throw std::invalid_argument("Precondition " #cond " failed")

// Postconditions check (results assertion)
#define Ensures(cond) if (!(cond)) throw std::runtime_error("Postcondition " #cond " failed")

// Runtime assertion. Why duplicate assert from cassert? Because we don't want to disable is in release mode and
// be more C++ and throw exception
#define Assert(cond) if (!(cond)) throw std::runtime_error("Assertion " #cond " failed")

// Additional macros for validating things like input from file - wrong input shouldn't be considered as assertion fail,
// because it is not programmer's fault ;)
struct ValidationException : public std::domain_error {
    explicit ValidationException(const std::string &msg) : std::domain_error{msg} { }
};

#define Validate(cond) if (!(cond)) throw ValidationException("Validation " #cond " failed")
#define ValidateMsg(cond, msg) if (!(cond)) throw ValidationException(msg)

#endif //RSA3D_ASSERTIONS_H
