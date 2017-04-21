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
	this->voxels = new VoxelList(dim, s, vdx);

	this->iAnalyze = 0;
	this->tmpSplit = 0;
	this->iMaxVoxels = 0;
	this->dMinVoxelSize = 0;
	this->missCounter = 0;
}

Surface::~Surface() {
	delete this->list;
	delete this->voxels;
}

void Surface::setParameters(int ia, int is, double dvs, int imv) {
	this->iAnalyze = ia;
	this->tmpSplit = is;
	this->iMaxVoxels = imv;
	this->dMinVoxelSize = dvs;
}

void Surface::setSeed(int s){
	this->seed = s;
}

void Surface::add(Shape *s) {
		s->no = this->shapes.size();
		this->shapes.insert(this->shapes.end(), s);
		this->list->add(s);
	}

bool Surface::check(Shape *s){
	bool bRet = true;
	std::unordered_set<Positioned *> * neighbours;
	neighbours = this->list->getNeighbours(s->getPosition());

	for(Positioned *shape: *neighbours) {
		if (((Shape *)shape)->overlap(this, s)){
			bRet = false;
			break;
		}
	}
	return bRet;
}


std::unordered_set<Positioned *> * Surface::getNeighbours(double* da) {
	return this->list->getNeighbours(da);
}

double Surface::distance2(double *a1, double *a2) {
	double v[this->dimension];
	for (int i = 0; i < this->dimension; i++)
		v[i] = a1[i] - a2[i];
	this->vector(v);
	double res = 0.0;
	for (int i = 0; i < this->dimension; i++)
		res += v[i] * v[i];
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

int Surface::analyzeVoxels() {
	std::cout << "[" << this->seed << "] " << this->voxels->length() << " voxels, " << this->shapes.size() << " shapes, factor = " << this->getFactor() << std::endl;
	int begin = this->voxels->length();
	int timestamp = this->shapes.size();
	for (int i = 0; i < this->voxels->length(); i++) {
		Voxel *v = this->voxels->get(i);
		if (this->voxels->analyzeVoxel(v, this->list, this, timestamp)) {
			this->voxels->remove(v);
			i--;
		}
	}
	std::cout << "[" << this->seed << "] " << this->voxels->length() << " voxels remained, factor = " << this->getFactor() << std::endl;
	return begin - this->voxels->length();
}

// analyzes all voxels inside a region around v
int Surface::analyzeRegion(Voxel *v){
	int begin = this->voxels->length();
	std::unordered_set<Positioned *> *region = this->voxels->getNeighbours(v);
	for(Positioned *v1: *region){
		if (this->voxels->analyzeVoxel((Voxel *)v1, this->list, this))
			this->voxels->remove((Voxel *)v1);
	}
	return begin - this->voxels->length();
}

bool Surface::doIteration(Shape *s, RND *rnd) {
	Voxel *v = this->voxels->getRandomVoxel(rnd);
	double da[this->dimension];
	s->translate(this->voxels->getRandomPosition(da, v, rnd));
	if (this->check(s)) {
		this->add(s);
		if(this->getFactor()>FACTOR_LIMIT){
			this->analyzeRegion(v);
		}else{
			this->voxels->remove(v);
		}
		this->missCounter = 0;
		return true;
	} else {
		v->miss();
		this->missCounter++;
		if (v->getMissCounter() % iAnalyze == 0) {
			if(this->voxels->analyzeVoxel(v, this->list, this) && this->getFactor()>FACTOR_LIMIT)
				this->analyzeRegion(v);
		}
		if (missCounter > tmpSplit) { // v.getMissCounter() % iSplit == 0){ //
			missCounter = 0;
			int v0 = this->voxels->length();
			bool b = voxels->splitVoxels(dMinVoxelSize, iMaxVoxels, this->list, this);
			int v1 = this->voxels->length();
			std::cout << "[" << this->seed << "] " << this->shapes.size() << " shapes, " << this->voxels->length() << " voxels, new voxel size " << voxels->getVoxelSize() << ", factor " << this->getFactor() << std::endl;
			if (b) {
				tmpSplit *=  ((double)v1 / v0);
				if (tmpSplit > 0.5 * std::numeric_limits<int>::max())
					tmpSplit = 0.5 * std::numeric_limits<int>::max();
				if(voxels->length()<1000 && tmpSplit>100)
					tmpSplit /= 10.0;
			} else {
				this->analyzeVoxels();
			}
		}
		return false;
	}
}

bool Surface::isSaturated() {
	return (this->voxels->length() == 0);
}

double Surface::getFactor() {
	return this->getArea() / this->voxels->getVoxelsSurface();
	}

std::vector<Positioned *> * Surface::getShapes() {
	return &(this->shapes);
}

