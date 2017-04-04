/*
 * Surface.cpp
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#include "Surface.h"
#include "NeighbourGrid.h"

#include <stdlib.h>

Surface::Surface(int dim, double s, double ndx, double vdx) {
	this->dimension = dim;
	this->size = s;
	this->list = new NeighbourGrid(dim, s, ndx);
//	this->voxels = new VoxelList(N, s, vdx);

	this->iAnalyze = 0;
	this->tmpSplit = 0;
	this->iMaxVoxels = 0;
	this->dMinVoxelSize = 0;
	this->missCounter = 0;
}

Surface::~Surface() {
	delete this->list;
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
	std::vector<Shape *> neighbours;
	if (this->list == NULL) {
		neighbours = this->shapes;
	} else {
		neighbours = this->list->getNeighbours(s->getPosition());
	}

	std::iterator = neighbours.begin();
	for(Shape * shape: neighbours) {
		if (shape->overlap(this, s))
			return false;
	}
	return true;
}


std::vector<Shape *> Surface::getNeighbours(double* da) {
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

/*
int Surface::analyzeVoxels() {
		std::cout << "[" << this->seed << "] " + this.voxels.length() + " voxels, " + this.shapes.size() + " shapes, factor = " + this.getFactor());
		int begin = this.voxels.length();
		int timestamp = this.shapes.size();
		for (int i = 0; i < this.voxels.length(); i++) {
			Voxel v = this.voxels.get(i);
			if (this.voxels.analyzeVoxel(v, this.list, this, timestamp)) {
				this.voxels.remove(v);
				i--;
			}
		}
		Globals.logger.info("[" + this.seed + "] " + this.voxels.length() + " voxels remained, factor = " + this.getFactor());
		return begin - this.voxels.length();
	}

	// analyzes all voxels inside a region around v
	private int analyzeRegion(Voxel v){
		int begin = this.voxels.length();
		Collection<Voxel> region = this.voxels.getNeighbours(v);
		for(Voxel v1: region){
			if (this.voxels.analyzeVoxel(v1, this.list, this))
				this.voxels.remove(v1);
			}
		return begin - this.voxels.length();
	}
*/


void Surface::setParameters(int ia, int is, double dvs, int imv) {
	this->iAnalyze = ia;
	this->tmpSplit = is;
	this->iMaxVoxels = imv;
	this->dMinVoxelSize = dvs;
}

/*
bool Surface::doIteration(Shape *s, RND *rnd) {
		Voxel v = this.voxels.getRandomVoxel(rnd);
		s.translate(this.voxels.getRandomPosition(v, rnd));
		if (this.check(s)) {
			this.add(s);
			if(this.getFactor()>FACTOR_LIMIT){
				this.analyzeRegion(v);
			}else{
				this.voxels.remove(v);
			}
			this.missCounter = 0;
			return true;
		} else {
			v.miss();
			this.missCounter++;
			if (v.getMissCounter() % iAnalyze == 0) {
//				this.voxels.analyzeVoxel(v, this.list, this, this.shapes.size());
				if(this.voxels.analyzeVoxel(v, this.list, this) && this.getFactor()>FACTOR_LIMIT)
					this.analyzeRegion(v);
			}
			if (missCounter > tmpSplit) { // v.getMissCounter() % iSplit == 0){ //
				missCounter = 0;
				int v0 = this.voxels.length();
				boolean b;
				if (s instanceof shapes.Disk || s instanceof shapes.NSphere) {
					b = voxels.splitVoxels(dMinVoxelSize, iMaxVoxels, this.list, this);
					Globals.logger.info("[" + this.seed + "] new voxel size " + voxels.getVoxelSize());
				} else {
					b = voxels.splitVoxels(dMinVoxelSize, iMaxVoxels, null, null);
					Globals.logger.info("[" + this.seed + "] new voxel size " + voxels.getVoxelSize());
					this.analyzeVoxels();
				}
				if (b) {
					if (this.getFactor() > FACTOR_LIMIT)
						s.shapeSpecificVoxelAnalysis(this, this.voxels);
					int v1 = this.voxels.length();
					if (v1 < v0)
						tmpSplit /= (2.0 * v0 / v1);
//					if (v1<10000)
//						tmpSplit = Math.min(tmpSplit,  iSplit);
				} else {
					this.analyzeVoxels();
				}
				if (tmpSplit < 0.4 * Integer.MAX_VALUE)
					tmpSplit *= 2.0;
			}
			return false;
		}
	}

	public boolean isSaturated() {
		return (this.voxels.length() == 0);
	}

	public double getFactor() {
		return this.getArea() / this.voxels.getVoxelsSurface();
	}

	public Shape[] getShapes() {
		return this.shapes.toArray(new Shape[0]);
	}

	synchronized public void drawShapes(Graphics g, double scale) {
		for (Shape s : this.shapes) {
			s.draw(g, scale);
		}
	}

	synchronized public void drawShapes(Graphics g, double scale, double[] ta) {
		for (Shape s : this.shapes) {
			Shape sc;
			try {
				sc = (Shape) s.clone();
				sc.translate(ta);
				sc.draw(g, scale);
			} catch (CloneNotSupportedException e) {
				e.printStackTrace();
			}
		}
	}
*/

