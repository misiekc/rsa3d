//
// Created by PKua on 21.07.18.
//

#ifndef RSA3D_ORDERCALCULATORS_H
#define RSA3D_ORDERCALCULATORS_H


#include <vector>
#include "../geometry/Vector.h"


class OrderParameters {
private:
    static double cosNSum(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes,
                          unsigned char cosExp);
    static double cosOfMaxAbs(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes);

public:
    static double nematicP1(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes);
    static double nematicP2(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes);
    static double nematicD2(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes);
    static double tetrahedral(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes);
    static double cubatic(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes);
    static double dodecahedral(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes);
};

namespace integer_legendre {
    inline double P1(double x) {
        return x;
    }

    inline double P2(double x) {
        return 0.5 * (3 * x * x - 1);
    }

    inline double P3(double x) {
        return 0.5 * (5 * x * x * x - 3 * x);
    }

    inline double P4(double x) {
        return 0.125 * (x * x * (35 * x * x - 30) + 3);
    }
}

#endif //RSA3D_ORDERCALCULATORS_H
