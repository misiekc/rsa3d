//
// Created by PKua on 21.07.18.
//

#ifndef RSA3D_ORDERCALCULATORS_H
#define RSA3D_ORDERCALCULATORS_H


#include <vector>
#include "../Vector.h"


class OrderParameters {
private:
    static double cosNSum(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes,
                          unsigned char cosExp);

public:
    static double nematic(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes);
    static double cubatic(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes);
    static double dodecahedral(const std::vector<Vector<3>> &firstAxes, const std::vector<Vector<3>> &secondAxes);
};


#endif //RSA3D_ORDERCALCULATORS_H
