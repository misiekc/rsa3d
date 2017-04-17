/*
 * Shape.cpp
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#include "Shape.h"
#include "Positioned.h"

Shape::Shape(const int dimension) : Positioned(dimension){
	this->no = 0;
	this->time = 0.0;
};

Shape::~Shape() {

}

void Shape::translate(double* v){
		for(int i=0; i<this->dimension; i++){
			this->position[i] += v[i];
		}
	}

