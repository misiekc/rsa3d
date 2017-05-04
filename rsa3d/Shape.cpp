/*
 * Shape.cpp
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#include "Shape.h"
#include "Positioned.h"

Shape::Shape(const unsigned short dimension) : Positioned(dimension){
	this->no = 0;
	this->time = 0.0;
}

Shape::Shape(const Shape & other) : Positioned(other){
	this->no = 0;
	this->time = 0.0;
}

Shape::~Shape() {

}

void Shape::translate(double* v){
	for(unsigned short i=0; i<this->dimension; i++){
		this->position[i] += v[i];
	}
}

std::string Shape::toString(){
	return "";
}

std::string Shape::toPovray(){
	return "";
}

void Shape::store(std::ostream &f){
}


