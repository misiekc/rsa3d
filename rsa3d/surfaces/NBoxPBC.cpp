/*
 * NBoxPBC.cpp
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#include "NBoxPBC.h"

NBoxPBC::NBoxPBC(double s, double ndx, double vdx) : Surface(s, ndx, vdx){
}

double NBoxPBC::getArea() const {
	return pow(this->size, RSA_SPATIAL_DIMENSION);
}

RSAVector NBoxPBC::getTranslation(double s, const RSAVector &p1, const RSAVector &p2) {
    RSAVector result;
	for(int i=0; i<RSA_SPATIAL_DIMENSION; i++){
		double d = 2*(p1[i] - p2[i]);
		if (d > s)
			result[i] = s;
		else if (d < -s)
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

void NBoxPBC::checkPosition(RSAVector &da) const {
	for(int i=0; i<RSA_SPATIAL_DIMENSION; i++) {
        if (da[i] < 0.0)
            da[i] += this->size;
        else if (da[i] >= this->size)
            da[i] -= this->size;
    }
}
