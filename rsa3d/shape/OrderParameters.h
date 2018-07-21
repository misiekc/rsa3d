//
// Created by PKua on 21.07.18.
//

#ifndef RSA3D_ORDERCALCULATORS_H
#define RSA3D_ORDERCALCULATORS_H


#include <vector>
#include "../Vector.h"

class OrderParameters {
public:
    struct OrderPair {
        double nematic;
        double full;
    };

    enum CosExp {
        _2 = 2,
        _4 = 4,
        _6 = 6
    };

    static OrderPair nematicAndFull(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes,
                                    CosExp cosExp);
    static double full(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes,
                       CosExp cosExp);
};


#endif //RSA3D_ORDERCALCULATORS_H
