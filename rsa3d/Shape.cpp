/*
 * Shape.cpp
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#include "Shape.h"

Shape::Shape(int dim) {
	this->dimension = dim;
	this->position = new double[dim];
}

Shape::~Shape() {
	delete this->position;
}


double* Shape::getPosition(){
	return this->position;

}

void Shape::translate(double* v){
		for(int i=0; i<this->dimension; i++){
			this->position[i] += v[i];
		}
	}

