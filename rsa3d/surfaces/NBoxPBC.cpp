/*
 * NBoxPBC.cpp
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#include "NBoxPBC.h"
#include "../utils/Assertions.h"

NBoxPBC::NBoxPBC(int dim, double s, double ndx, double vdx) : Surface(dim, s, ndx, vdx){
    ValidateMsg(s > 4*RSAShape::getCircumsphereRadius(),
                "Packing linear size is <= 4 circumscribed spheres - boundary conditions will break. ");
}

RSAVector NBoxPBC::getTranslation(double s, const RSAVector &p1, const RSAVector &p2) {
    RSAVector result;
    RSAVector delta = p1 - p2;
	for(int i=0; i<RSA_SPATIAL_DIMENSION; i++){
		if (delta[i] > s/2)
			result[i] = s;
		else if (delta[i] < -s/2)
			result[i] = -s;
	}
	return result;
}

RSAVector NBoxPBC::getTranslation(const RSAVector &p1, const RSAVector &p2) const {
    return NBoxPBC::getTranslation(this->size, p1, p2);
}

RSAVector NBoxPBC::vector(const RSAVector &v) const {
	return this->vectorPeriodicBC(v);
}

RSAVector NBoxPBC::checkPosition(const RSAVector &da) const {
    RSAVector result = da;
	for(int i=0; i<RSA_SPATIAL_DIMENSION; i++) {
        if (da[i] < 0.0)
            result[i] += this->size;
        else if (da[i] >= this->size)
            result[i] -= this->size;
    }
    return result;
}
