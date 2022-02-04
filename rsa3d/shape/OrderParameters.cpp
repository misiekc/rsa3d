//
// Created by PKua on 21.07.18.
//

#include "OrderParameters.h"
#include "../utils/Utils.h"

inline double OrderParameters::cosNSum(const std::vector<Vector<3>> &firstAxes,
                                       const std::vector<Vector<3>> &secondAxes, unsigned char cosExp) {
    double cosNsum = 0;
    for (const auto &a1 : firstAxes)
        for (const auto &a2 : secondAxes)
            cosNsum += std::pow(a1 * a2, cosExp);

    return cosNsum;
}

inline double OrderParameters::cosOfMaxAbs(const std::vector<Vector<3>> &firstAxes,
                                           const std::vector<Vector<3>> &secondAxes) {
    double maxAbsCos = 0;
    double maxCos = 0;
    for (const auto &a1 : firstAxes) {
        for (const auto &a2 : secondAxes) {
            double cos = a1 * a2;
            double absCos = std::abs(cos);
            if (absCos > maxAbsCos) {
                maxAbsCos = absCos;
                maxCos = cos;
            }
        }
    }

    return maxCos;
}

double OrderParameters::nematicP1(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes) {
    return integer_legendre::P1(cosOfMaxAbs(firstAxes, secondAxes));
}

double OrderParameters::nematicP2(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes) {
    return integer_legendre::P2(cosOfMaxAbs(firstAxes, secondAxes));
}

double OrderParameters::tetrahedral(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes) {
    return 9./32*cosNSum(firstAxes, secondAxes, 3);
}

double OrderParameters::cubatic(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes) {
    return 1./6*(5*cosNSum(firstAxes, secondAxes, 4) - 9);
}

double OrderParameters::dodecahedral(const std::vector<Vector<3>> &firstAxes,
                                     const std::vector<Vector<3>> &secondAxes) {
    return 25./192*(7*cosNSum(firstAxes, secondAxes, 6) - 36);
}
