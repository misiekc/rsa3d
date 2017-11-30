/*
 * Shape.cpp
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#include <iostream>

template <unsigned short DIMENSION>
Shape<DIMENSION>::Shape() : Positioned<DIMENSION>(){
	this->no = 0;
	this->time = 0.0;
}

template <unsigned short DIMENSION>
Shape<DIMENSION>::Shape(const Shape<DIMENSION> & other) : Positioned<DIMENSION>(other){
	this->no = 0;
	this->time = 0.0;
}

template <unsigned short DIMENSION>
Shape<DIMENSION>::~Shape() {

}

template <unsigned short DIMENSION>
void Shape<DIMENSION>::translate(double* v){
	for(unsigned short i=0; i<DIMENSION; i++){
		this->position[i] += v[i];
	}
}

template <unsigned short DIMENSION>
double Shape<DIMENSION>::minDistance(Shape *s){
	return 0.0;
}

template <unsigned short DIMENSION>
std::string Shape<DIMENSION>::toString(){
	return "";
}

template <unsigned short DIMENSION>
std::string Shape<DIMENSION>::toPovray() const{
	return "";
}

template <unsigned short DIMENSION>
void Shape<DIMENSION>::store(std::ostream &f) const{
	unsigned short dim = DIMENSION;
	f.write((char *)(&dim), sizeof(unsigned char));
	f.write((char *)(&this->no), sizeof(int));
	f.write((char *)(&this->time), sizeof(double));
	f.write((char *)(this->position), DIMENSION*sizeof(double));
}

template <unsigned short DIMENSION>
void Shape<DIMENSION>::restore(std::istream &f){
	unsigned char d;
	f.read((char *)(&d), sizeof(unsigned char));

	if (f.gcount()==0){ // end of file
		return;
	}
	if (d!=DIMENSION){
		std::cout << "[ERROR] cannot restore: incompatible dimensions: read " << f.gcount() << " bytes." << std::endl;
		return;
	}
	f.read((char *)(&this->no), sizeof(int));
	f.read((char *)(&this->time), sizeof(double));
	f.read((char *)(this->position), d*sizeof(double));
}


