//
// Created by pkua on 24.11.2021.
//

#include <catch.hpp>

#include "../../../../../rsa3d/shape/shapes/polygon/RoundedPolygon.h"

namespace {
    double get_segment_length(const std::vector<double> &vertexR, const std::vector<double> &vertexTheta,
                              std::size_t s1, std::size_t s2)
    {
        Vector<2> vertex1 = {{vertexR[s1] * std::cos(vertexTheta[s1]), vertexR[s1] * std::sin(vertexTheta[s1])}};
        Vector<2> vertex2 = {{vertexR[s2] * std::cos(vertexTheta[s2]), vertexR[s2] * std::sin(vertexTheta[s2])}};
        return (vertex2 - vertex1).norm();
    }
}

TEST_CASE("RoundedPolygon: area normalization") {
    SECTION("perfect square") {
        RoundedPolygon::initClass("1 4 xy 0.5 0.5 -0.5 0.5 -0.5 -0.5 0.5 -0.5 4 0 1 2 3");

        const auto &vertexR = Polygon::getVertexR();
        const auto &vertexTheta = Polygon::getVertexTheta();

        const double expectedLength = 1 / std::sqrt(5 + M_PI);
        CHECK(vertexR.size() == 4);
        CHECK(vertexTheta.size() == 4);
        CHECK(get_segment_length(vertexR, vertexTheta, 0, 1) == Approx(expectedLength));
        CHECK(get_segment_length(vertexR, vertexTheta, 1, 2) == Approx(expectedLength));
        CHECK(get_segment_length(vertexR, vertexTheta, 2, 3) == Approx(expectedLength));
        CHECK(get_segment_length(vertexR, vertexTheta, 3, 0) == Approx(expectedLength));
    }

    SECTION("obtuse triangle") {
        RoundedPolygon::initClass("0.405684 3 xy "
                                  "0.1784267738602887 -0.3683330410366601 "
                                  "-0.1027397233759038 -0.03388184324595755 "
                                  "-0.07568705048438495 0.4022148842826159 "
                                  "3 0 1 2");

        const auto &vertexR = Polygon::getVertexR();
        const auto &vertexTheta = Polygon::getVertexTheta();

        CHECK(vertexR.size() == 3);
        CHECK(vertexTheta.size() == 3);
        CHECK(get_segment_length(vertexR, vertexTheta, 0, 1) == Approx(0.3882451851870135));
        CHECK(get_segment_length(vertexR, vertexTheta, 1, 2) == Approx(0.3882451851870145));
        CHECK(get_segment_length(vertexR, vertexTheta, 2, 0) == Approx(0.7209532095891463));
    }
}