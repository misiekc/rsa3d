//
// Created by PKua on 21.07.18.
//

#include "OrderParameters.h"
#include "../Utils.h"

OrderParameters::OrderPair
OrderParameters::nematicAndFull(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes,
                                CosExp cosExp) {
    double cosNsum = 0;
    double maxAbsCos = 0;
    for (const auto &a1 : firstAxes) {
        for (const auto &a2 : secondAxes) {
            double absCos = std::abs(a1 * a2);
            if (absCos > maxAbsCos)
                maxAbsCos = absCos;
            cosNsum += std::pow(absCos, cosExp);
        }
    }

    double nematicOrder = P2(maxAbsCos);
    return {nematicOrder, cosNsum};
}

double OrderParameters::full(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes,
                             CosExp cosExp) {
    double cosNsum = 0;
    for (const auto &a1 : firstAxes)
        for (const auto &a2 : secondAxes)
            cosNsum += std::pow(a1 * a2, cosExp);

    return cosNsum;
}