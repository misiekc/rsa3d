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
	this->list = new NeighbourGrid(dim, s, ndx);
}

Surface::~Surface() {
	delete this->list;
}


void Surface::add(Shape *s) {
		this->list->add(s);
	}

Shape* Surface::check(Shape *s){
	std::unordered_set<Positioned *> neighbours;
	this->list->getNeighbours(&neighbours, s->getPosition());

	for(Positioned *shape: neighbours) {
		if (((Shape *)shape)->overlap(this, s)){
			return (Shape *)shape;
		}
	}
	return NULL;
}

Positioned * Surface::getClosestNeighbour(double *da){
	return this->list->getClosestNeighbour(da, this);
}

void Surface::getNeighbours(std::unordered_set<Positioned *> *result, double* da) {
	this->list->getNeighbours(result, da);
}

NeighbourGrid * Surface::getNeighbourGrid(){
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
