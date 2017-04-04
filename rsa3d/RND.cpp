/*
 * RND.cpp
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#include "RND.h"

RND::RND(int seed) {
	this->mt(seed);
	this->distribution(0.0, 1.0);
}

RND::RND() {
	this->mt();
	this->distribution(0.0, 1.0);
}

RND::~RND() {

}

double RND::nextValue(){
	return this->distribution(mt);
}

