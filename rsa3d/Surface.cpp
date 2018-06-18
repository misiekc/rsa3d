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
	this->list = new NeighbourGrid<const RSAShape>(dim, s, ndx);
}

Surface::~Surface() {
	delete this->list;
}

void Surface::clear(){
	this->list->clear();
}


void Surface::add(const RSAShape *s) {
	this->list->add(s, s->getPosition());
}

const RSAShape* Surface::check(const RSAShape *s){
	std::vector<const RSAShape *> neighbours;
	this->list->getNeighbours(&neighbours, s->getPosition());

	for(const RSAPositioned *positioned: neighbours) {
	    auto shape = dynamic_cast<const RSAShape*>(positioned);
		if (shape->overlap(this, s))
			return shape;
	}
	return nullptr;
}


const RSAShape *Surface::getClosestNeighbour(double *da, const std::vector<const RSAShape*> &neighbours){
    double d, dmin = std::numeric_limits<double>::max();
    const RSAShape *pmin = nullptr;
    for(const RSAShape *p : neighbours){
        d = this->distance2(da, p->getPosition());
        if (d<dmin){
            pmin = p;
            dmin = d;
        }
    }
    return pmin;
}


const RSAShape *Surface::getClosestNeighbour(double *da) {
    std::vector<const RSAShape*> result;
    this->list->getNeighbours(&result, da);
    return getClosestNeighbour(da, result);
}


void Surface::getNeighbours(std::vector<const RSAShape*> *result, double* da) {
	this->list->getNeighbours(result, da);
}

NeighbourGrid<const RSAShape> *Surface::getNeighbourGrid(){
	return this->list;
}

void Surface::checkPosition(double *da){
	// do nothing
}

double Surface::distance2(const double *a1, const double *a2) {
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
