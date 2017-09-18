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

NBoxFBC::~NBoxFBC() {
	// TODO Auto-generated destructor stub
}

double NBoxFBC::getArea(){
	return pow(this->size, this->dimension);
}

double * NBoxFBC::getTranslation(double *result, int dim, double s, double *p1, double *p2) {
	for(int i=0; i<dim; i++)
		result[i] = 0.0;
	return result;
}

double * NBoxFBC::getTranslation(double *result, double *p1, double *p2) {
		return NBoxFBC::getTranslation(result, this->dimension, this->size, p1, p2);
	}

void NBoxFBC::vector(double* v) {
	this->vectorFreeBC(v);
}
