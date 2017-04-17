/*
 * PackingGenerator.cpp
 *
 *  Created on: 16.04.2017
 *      Author: ciesla
 */

#include "PackingGenerator.h"
#include "RND.h"
#include "surfaces/NBoxPBC.h"
#include "ShapeFactory.h"
#include <iostream>

PackingGenerator::PackingGenerator(int s) {
	this->params = new Parameters();
	this->seed = s;

}

PackingGenerator::~PackingGenerator() {
	delete this->params;
	for(Shape *s : this->packing)
		delete s;
}

void PackingGenerator::createPacking(){
	std::cout << "[" << this->seed << "] started" << std::endl;
	int missCount = 0;
	RND rnd(this->seed);
	ShapeFactory::initShapeClass(this->params->particleType, this->params->particleAttributes);
	Shape *s = ShapeFactory::createShape(&rnd);
	Surface *surface = new NBoxPBC(this->params->dimension, this->params->surfaceSize, s->getNeighbourListCellSize(), s->getVoxelSize());
	surface->setSeed(this->seed);
	double dt = s->getVolume() / surface->getArea();
	delete s;
	int l = 0;
	double t = 0;
	surface->setParameters(params->analyze, params->split, params->minDx, params->maxVoxels);
	while (!surface->isSaturated() && t<params->maxTime && missCount<params->maxTriesWithoutSuccess) {
		t += surface->getFactor() * dt;
		s = ShapeFactory::createShape(&rnd);
		if (surface->doIteration(s, &rnd)) {
			l++;
			s->no = l;
			s->time = t;
			this->packing.push_back(s);
			// double[] da = s.getCoordinates();
			if (t>0.1*params->maxTime)
				std::cout << "[" << this->seed << "]" << "\t" << t << "\t" << surface->getFactor()
				<< "\t" << l << "\t" << surface->voxels->length()
				<< "\t" << missCount << std::endl;
			missCount = 0;
		}else{
			missCount++;
		}
	}
	delete surface;
	std::cout << "[" << seed << "] finished after generating " << l << " shapes" << std::endl;
}

void PackingGenerator::run(){
	this->createPacking();
	for(Shape *s : this->packing){
		double *da = s->getPosition();
		std::cout << s->no << "\t";
		for(int i=0; i<this->params->dimension; i++)
			std::cout << da[i] << "\t";
		std::cout << std::endl;
	}
}
