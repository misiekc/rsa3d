/*
 * Surface.cpp
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#include "Surface.h"
#include <iostream>

Surface::Surface(int dim, double s, double ndx, double vdx, std::unique_ptr<RSABoundaryConditions> bc) {
    ValidateMsg(s > 2*RSAShape::getNeighbourListCellSize(),
                "Packing linear size is <= 2 neighbour list cell size - boundary conditions will break. ");

    this->bc = std::move(bc);
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

double Surface::getArea() const {
    return pow(this->size, this->dimension);
}

double Surface::getRealArea() const {
    return pow(this->size, this->dimension);
}

void Surface::add(const RSAShape *s) {
	this->list->add(s, s->getPosition());
}

const RSAShape* Surface::check(const RSAShape *s){
	std::vector<const RSAShape *> neighbours;
	this->list->getNeighbours(&neighbours, s->getPosition());
	return s->overlap(this->bc.get(), &neighbours);
/*
	for(const RSAShape *shape: neighbours) {
		if (shape->overlap(this, s))
			return shape;
	}
	return nullptr;
*/
}


const RSAShape *Surface::getClosestNeighbour(const RSAVector &da, const std::vector<const RSAShape*> &neighbours){
    double d, dmin = std::numeric_limits<double>::max();
    const RSAShape *pmin = nullptr;
    for(const RSAShape *p : neighbours){
        d = this->bc->distance2(da, p->getPosition());
        if (d<dmin){
            pmin = p;
            dmin = d;
        }
    }
    return pmin;
}

const RSAShape *Surface::getClosestNeighbour(const RSAVector &da) {
    std::vector<const RSAShape*> result;
    this->list->getNeighbours(&result, da);
    return getClosestNeighbour(da, result);
}


void Surface::getNeighbours(std::vector<const RSAShape*> *result, const RSAVector &da) {
	this->list->getNeighbours(result, da);
}

NeighbourGrid<const RSAShape> *Surface::getNeighbourGrid(){
	return this->list;
}

RSAVector Surface::checkPosition(const RSAVector &da) const {
	return da;	// do nothing
}


RSABoundaryConditions *Surface::getBC() {
    return this->bc.get();
}
