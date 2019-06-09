/*
 * ArrayFunction.cpp
 *
 *  Created on: Jun 7, 2019
 *      Author: ciesla
 */

#include "ArrayFunction.h"



ArrayFunction::ArrayFunction(double **data, unsigned int length) {
	this->data = new double*[length];
	for(unsigned int i=0; i<length; i++){
		this->data[i] = new double[2];
		this->data[i][0] = data[i][0];
		this->data[i][1] = data[i][1];
	}
	this->length = length;
}

ArrayFunction::~ArrayFunction() {
	for(unsigned int i=0; i<length; i++)
		delete[] this->data[i];
	delete[] this->data;
}

double ArrayFunction::get(double x){
	double x1, x2, y1, y2;
	int i1, i2;
	if (x<this->data[0][0]){
		i1 = 0; i2 =1;
	}else if (x>this->data[this->length-1][0]){
		i1 = this->length-2;
		i2 = this->length-1;
	}else{
		int i;

		i1 = 0;
		i2 = this->length;
		while(i2-i1>1){
			i = (i1+i2)/2;
			if (this->data[i][0] < x)
				i1 = i;
			else
				i2 = i;
		}
		if (i1==i2){
			i2++;
		}
	}
	x1 = this->data[i1][0]; x2 = this->data[i2][0];
	y1 = this->data[i1][1]; y2 = this->data[i2][1];
	return y1 + (x-x1)*(y2-y1)/(x2-x1);
}

