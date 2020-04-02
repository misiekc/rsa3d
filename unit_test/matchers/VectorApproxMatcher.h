//
// Created by Piotr Kubala on 27/03/2020.
//

#ifndef RSA3D_VECTORAPPROXMATCHER_H
#define RSA3D_VECTORAPPROXMATCHER_H

#include <catch.hpp>
#include <sstream>
#include <iterator>

#include "../../rsa3d/geometry/Vector.h"

template<std::size_t DIM>
class VectorApproxEqualMatcher : public Catch::MatcherBase<Vector<DIM>> {
private:
    Vector<DIM> expected;
    double epsilon;

public:
    VectorApproxEqualMatcher(Vector<DIM> expected, double epsilon) : expected(std::move(expected)), epsilon(epsilon) { }

    bool match(const Vector<2> &actual) const override {
        for (std::size_t i{}; i < DIM; i++)
            if (this->expected[i] != Approx(actual[i]).epsilon(this->epsilon))
                return false;
        return true;
    }

    [[nodiscard]] std::string describe() const override {
        std::ostringstream ss;
        ss << "is, within " << epsilon << " tolerance threshold, equal to " << this->expected;
        return ss.str();
    }
};

template<std::size_t DIM>
inline VectorApproxEqualMatcher<DIM> IsApproxEqual(const Vector<DIM> &expected, double epsilon) {
    return VectorApproxEqualMatcher(expected, epsilon);
}


#endif //RSA3D_VECTORAPPROXMATCHER_H
