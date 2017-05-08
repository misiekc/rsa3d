/*
 * NBoxPBC.cpp
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#include "NBoxPBC.h"

NBoxPBC::NBoxPBC(int dim, double s, double ndx, double vdx) : Surface(dim, s, ndx, vdx){
}

NBoxPBC::~NBoxPBC() {
}

double NBoxPBC::getArea(){
	return pow(this->size, this->dimension);
}

double * NBoxPBC::getTranslation(double *result, int dim, double s, double *p1, double *p2) {
	double d;
	for(int i=0; i<dim; i++){
		d = 2*(p1[i] - p2[i]);
		if (d > s)
			result[i] = s;
		else if (d < -s)
			result[i] = -s;
		else
			result[i] = 0.0;
		}
	return result;
	}

double * NBoxPBC::getTranslation(double *result, double *p1, double *p2) {
		return NBoxPBC::getTranslation(result, this->dimension, this->size, p1, p2);
	}

void NBoxPBC::vector(double* v) {
	this->vectorPeriodicBC(v);
}

double * NBoxPBC::getRandomPosition(double *result, RND *rnd) {
	for(int i=0; i<this->dimension; i++)
		result[i] = rnd->nextValue()*this->size;
	return result;
}

