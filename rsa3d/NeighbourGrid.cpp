/*
 * NeighbourGrid.cpp
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#include "NeighbourGrid.h"
#include <cmath>
#include <algorithm>

#include "Utils.h"

/**
 *
 * @param dim dimension
 * @param size linear size of the structure
 * @param dx linear size of one cell
 */
NeighbourGrid::NeighbourGrid(int dim, double size, double dx) {
	this->dimension = dim;
	this->linearSize = size;
	this->n = (int)(size/dx);
	this->dx = size/this->n;
	int length = (int) round(pow(this->n, dim));
	this->lists.reserve(length);
	for(int i=0; i<length; i++){
		this->lists.push_back(new std::vector<Positioned *>);
	}
	this->neighbours.clear();

	this->in = new int[dim];
}


NeighbourGrid::~NeighbourGrid() {
	for(unsigned int i=0; i<this->lists.size(); i++){
		delete this->lists[i];
	}

	delete[] this->in;
}

void NeighbourGrid::add(Positioned* s){
	double* da = s->getPosition();
	int i = position2i(da, this->dimension, this->linearSize, this->dx, this->n);
	this->lists[i]->push_back(s);
}

/*

	protected static <T> boolean compareLists(ArrayList<T> l1, ArrayList<T> l2){
		if(l1.size()!=l2.size())
			return false;
		for(T t: l1){
			if (!l2.contains(t))
				return false;
		}
		return true;
	}
*/


void NeighbourGrid::remove(Positioned* s){
	double* da = s->getPosition();
	int i = position2i(da, this->dimension, this->linearSize, this->dx, this->n);
	std::vector<Positioned *>::iterator it;
	if ( (it = std::find(this->lists[i]->begin(), this->lists[i]->end(), s)) != this->lists[i]->end())
		this->lists[i]->erase(it);
}

std::unordered_set<Positioned*> * NeighbourGrid::getNeighbours(double* da, int radius){
	this->neighbours.clear();
	std::vector<Positioned *> *vTmp;

	for(int i=0; i<this->dimension; i++){
		this->in[i] = 0;
	}

	int coords[this->dimension];

	coordinates(coords, da, this->dimension, this->linearSize, this->dx, this->n);
	do{
		int i = neighbour2i(coords, this->in, this->dimension, radius, this->n);
		vTmp = (this->lists[i]);
		this->neighbours.insert(vTmp->begin(), vTmp->end());
	}while(increment(in, this->dimension, 2*radius));
	return &this->neighbours;
}

void NeighbourGrid::clear(){
	for(unsigned int i=0; i<this->lists.size(); i++){
		this->lists[i]->clear();
	}
}

std::unordered_set<Positioned*> * NeighbourGrid::getNeighbours(double* da){
	return this->getNeighbours(da, 1);
}

