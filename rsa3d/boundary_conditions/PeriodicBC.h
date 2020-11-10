/*
 * NBoxPBC.h
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#ifndef BOUNDARY_CONDITIONS_PERIODICBC_H_
#define BOUNDARY_CONDITIONS_PERIODICBC_H_

#include "../BoundaryConditions.h"
#include "../utils/Assertions.h"

template<unsigned short DIMENSION>
class PeriodicBC : public RSABoundaryConditions {
private:
    using vector = Vector<DIMENSION>;

    double surfaceSize{};

    static vector getTranslation(double s, const vector &p1, const vector &p2) {
        vector result;
        vector delta = p1 - p2;
        for(int i=0; i<DIMENSION; i++){
            if (delta[i] > s/2)
                result[i] = s;
            else if (delta[i] < -s/2)
                result[i] = -s;
        }
        return result;
    }

public:
    explicit PeriodicBC(double surfaceSize) : surfaceSize{surfaceSize} {
        Expects(surfaceSize > 0);
    }

    vector getTranslation(const vector &p1, const vector &p2) const override {
        return PeriodicBC::getTranslation(this->surfaceSize, p1, p2);
    }

    [[nodiscard]] double distance2(const vector &p1, const vector &p2) const override {
        vector distance = p2 - p1;
        for (int i = 0; i < DIMENSION; i++) {
            if (distance[i] > this->surfaceSize / 2.0)
                distance[i] -= this->surfaceSize;
            else if (distance[i] < -this->surfaceSize / 2.0)
                distance[i] += this->surfaceSize;
        }
        return distance.norm2();
    }
};

using RSAPeriodicBC = PeriodicBC<RSA_SPATIAL_DIMENSION>;

#endif /* BOUNDARY_CONDITIONS_PERIODICBC_H_ */
