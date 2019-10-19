/*
 * NBoxFBC.cpp
 *
 *  Created on: 13.07.2017
 *      Author: ciesla
 */

#include "NBoxFBC.h"

NBoxFBC::NBoxFBC(int dim, double s, double ndx, double vdx) : Surface(dim, s, ndx, vdx) {
	// TODO Auto-generated constructor stub
}

RSAVector NBoxFBC::getTranslation(double s, const RSAVector &p1, const RSAVector &p2) {
	return {};
}

RSAVector NBoxFBC::getTranslation(const RSAVector &p1, const RSAVector &p2) const {
	return NBoxFBC::getTranslation(this->size, p1, p2); }

RSAVector NBoxFBC::vector(const RSAVector &v) const {
	return this->vectorFreeBC(v);
}
