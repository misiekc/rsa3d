/*
 * Surface.cpp
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#include "Surface.h"
#include <iostream>

Surface::Surface(int dim, double s, double ndx, double vdx) : BoundaryConditions() {
	this->dimension = dim;
	this->size = s;
	this->list = new NeighbourGrid<Shape<RSA_DIMENSION>>(dim, s, ndx);
}

Surface::~Surface() {
	delete this->list;
}


void Surface::add(Shape<RSA_DIMENSION> *s) {
		this->list->add(s, s->getPosition());
	}

Shape<RSA_DIMENSION>* Surface::check(Shape<RSA_DIMENSION> *s){
	std::vector<Shape<RSA_DIMENSION> *> neighbours;
	this->list->getNeighbours(&neighbours, s->getPosition());

	for(Positioned<RSA_DIMENSION> *shape: neighbours) {
		if (((Shape<RSA_DIMENSION> *)shape)->overlap(this, s)){
			return (Shape<RSA_DIMENSION> *)shape;
		}
	}
	return NULL;
}


Shape<RSA_DIMENSION> * Surface::getClosestNeighbour(double *da, std::vector<Shape<RSA_DIMENSION> *> *neighbours){

	std::vector<Shape<RSA_DIMENSION> *> result;
	if(neighbours==NULL){
		this->list->getNeighbours(&result, da);
		neighbours = &result;
	}
	double d, dmin = std::numeric_limits<double>::max();
	Shape<RSA_DIMENSION> *pmin = NULL;
	for(Shape<RSA_DIMENSION> *p : *neighbours){
		d = this->distance2(da, p->getPosition());
		if (d<dmin){
			pmin = p;
			dmin = d;
		}
	}
	return pmin;
}



void Surface::getNeighbours(std::vector<Shape<RSA_DIMENSION> *> *result, double* da) {
	this->list->getNeighbours(result, da);
}

NeighbourGrid<Shape<RSA_DIMENSION>> * Surface::getNeighbourGrid(){
	return this->list;
}

void Surface::checkPosition(double *da){
	// do nothing
}

bool Surface::isInside(double *da){
	for (int i = 0; i < this->dimension; i++)
		if (da[i] < 0 || da[i] > this->size)
			return false;
	return true;
}

double Surface::distance2(double *a1, double *a2) {
	double *v = new double[this->dimension];
	for (int i = 0; i < this->dimension; i++)
		v[i] = a1[i] - a2[i];
	this->vector(v);
	double res = 0.0;
	for (int i = 0; i < this->dimension; i++)
		res += v[i] * v[i];
	delete[] v;
	return res;
}

void Surface::vectorFreeBC(double* v) {
	// do nothing
}

void Surface::vectorPeriodicBC(double* v) {
	for (int i = 0; i < this->dimension; i++) {
		if (v[i] > this->size / 2.0)
			v[i] -= this->size;
		else if (v[i] < -this->size / 2.0)
			v[i] += this->size;
	}
}
