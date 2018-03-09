/*
 * PackingGenerator.cpp
 *
 *  Created on: 16.04.2017
 *      Author: ciesla
 */

#include "RND.h"
#include "surfaces/NBoxPBC.h"
#include "surfaces/NBoxFBC.h"
#include "ShapeFactory.h"
#include <iostream>
#include <fstream>

#include "shapes/Ellipse.h"


#ifdef _OPENMP
#include <omp.h>
#endif

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::FACTOR_LIMIT = 5.0;

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::PackingGenerator(int s, Parameters *p) {
	this->params = p;
	this->seed = s;
	this->voxels = NULL;
	this->surface = NULL;

}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::~PackingGenerator() {
	for(Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s : this->packing)
		delete s;
	if (this->voxels!=NULL)
		delete this->voxels;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::isSaturated() {
	return (this->voxels->length() == 0);
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getFactor() {
	return this->surface->getArea() / this->voxels->getVoxelsSurface();
}

#ifdef _OPENMP
template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
int PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::analyzeVoxels() {

	std::cout << "[" << this->seed << " PackingGenerator::analyzeVoxels] " << this->voxels->length() << " voxels, " << this->packing.size() << " shapes, factor = " << this->getFactor() << std::endl;

	int begin = this->voxels->length();
	int timestamp = this->packing.size();
	std::vector< Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>* > toRemove;

	#pragma omp parallel for
	for (int i = 0; i < this->voxels->length(); i++) {
		Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v = this->voxels->get(i);
		if (this->voxels->analyzeVoxel(v, this->surface->getNeighbourGrid(), this->surface, timestamp)) {
			#pragma omp critical
			toRemove.push_back(v);
		}
	}

	for(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v: toRemove){
		this->voxels->remove(v);
	}

	std::cout << "[" << this->seed << " Surface::analyzeVoxels] " << this->voxels->length() << " voxels remained, factor = " << this->getFactor() << std::endl;

	return begin - this->voxels->length();
}

#else

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
int PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::analyzeVoxels() {

	std::cout << "[" << this->seed << " PackingGenerator::analyzeVoxels] " << this->voxels->length() << " voxels, " << this->packing.size() << " shapes, factor = " << this->getFactor() << std::endl;

	int begin = this->voxels->length();
	int timestamp = this->packing.size();

	for (int i = 0; i < this->voxels->length(); i++) {
		Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v = this->voxels->get(i);
		if (this->voxels->analyzeVoxel(v, this->surface->getNeighbourGrid(), this->surface, timestamp)) {
			this->voxels->remove(v);
			i--;
		}
	}

	std::cout << "[" << this->seed << " Surface::analyzeVoxels] " << this->voxels->length() << " voxels remained, factor = " << this->getFactor() << std::endl;

	return begin - this->voxels->length();
}

#endif

// analyzes all voxels inside a region around v
template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
int PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::analyzeRegion(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v){
	int begin = this->voxels->length();
	std::vector<Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> region;
	this->voxels->getNeighbours(&region, v);
	for(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v1: region){
		if (this->voxels->analyzeVoxel(v1, this->surface->getNeighbourGrid(), this->surface))
			this->voxels->remove(v1);
	}
	return begin - this->voxels->length();
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::modifiedRSA(Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s, Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v){
	double da[SPATIAL_DIMENSION];

	Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *sn = (Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *)this->surface->getClosestNeighbour(s->getPosition(), NULL);
	if (sn==NULL)
		sn = (Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *)this->surface->getClosestNeighbour(s->getPosition(), &(this->packing));
	if (sn!=NULL){
		double *spos = s->getPosition();
		double *snpos = sn->getPosition();

		double d = sqrt(this->surface->distance2(spos, snpos));
		if (d < this->params->thresholdDistance){

			for(ushort i=0; i<SPATIAL_DIMENSION; i++){
				da[i] = (snpos[i] - spos[i]);
			}
			this->surface->vector(da);
			double mindist = s->minDistance(sn);
			for(ushort i=0; i<SPATIAL_DIMENSION; i++){
				da[i] = da[i]/d;
				spos[i] += (d-mindist)*da[i];
			}
			this->surface->checkPosition(spos);
			v = this->voxels->getVoxel(spos, s->getOrientation());
			if (v==NULL){
				std::cout << "Problem: PackingGenerator - voxel not found: " <<
				" (" << spos[0] << ", " << spos[1] << ")" << std::endl;
				exit(0);
			}
		}
	}
}

#ifdef _OPENMP
template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::createPacking(){

	std::cout << "[" << this->seed << " PackingGenerator::createPacking] started" << std::endl;

	int missCounter = 0;
	RND rnd(this->seed);
	Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s = ShapeFactory::createShape(&rnd);

	this->voxels = new VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>(this->params->surfaceSize, s->getVoxelSize());

	double gridSize = s->getNeighbourListCellSize();
	if (gridSize < this->params->thresholdDistance)
		gridSize = this->params->thresholdDistance;


	if(params->boundaryConditions.compare("free")==0)
		this->surface = new NBoxFBC(SPATIAL_DIMENSION, this->params->surfaceSize, s->getNeighbourListCellSize(), s->getVoxelSize());
	else
		this->surface = new NBoxPBC(SPATIAL_DIMENSION, this->params->surfaceSize, s->getNeighbourListCellSize(), s->getVoxelSize());

	double dt = s->getVolume() / this->surface->getArea();
	delete s;
	int l = 0;
	double t = 0;
	int tmpSplit = this->params->split;
//	int snapshotCounter = 0;
	RND **aRND = NULL;
	int aRNDSize = 0;
	Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> **sOverlapped = new Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *[tmpSplit];
	Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> **sVirtual = new Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *[tmpSplit];
	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> **aVoxels = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *[tmpSplit];

	while (!this->isSaturated() && t<params->maxTime && missCounter<params->maxTriesWithoutSuccess) {

		#pragma omp parallel
		{

			#pragma omp master
			{
//				std::cout << "[master] init aRND[" << omp_get_num_threads() << "]" << std::endl << std::flush;

				if (aRND==NULL){ // if there are no random number generators
					aRND = new RND*[omp_get_num_threads()];
					for(int i=0; i<omp_get_num_threads(); i++){
						aRND[i] = new RND((int)(1000*rnd.nextValue()));
					}
				}else if(omp_get_num_threads() > aRNDSize){ // if the number of allocated random number generators is smaller than the number of threads
					for(int i=aRNDSize; i<omp_get_num_threads(); i++){
						aRND[i] = new RND((int)(1000*rnd.nextValue()));
					}
				}
				aRNDSize = omp_get_num_threads();
			}

			#pragma omp barrier

//			#pragma omp critical
//			std::cout << "[" << omp_get_thread_num() << "] looping" << std::endl << std::flush;

			#pragma omp for
			for(int i = 0; i<tmpSplit; i++){
				int tid = omp_get_thread_num();

				sVirtual[i] = ShapeFactory::createShape(aRND[tid]);
				double pos[SPATIAL_DIMENSION];
				double angle[ANGULAR_DIMENSION];
				do{
					aVoxels[i] = this->voxels->getRandomVoxel(aRND[tid]);
					this->voxels->getRandomPositionAndOrientation(pos, angle, aVoxels[i], aRND[tid]);
				}while(!this->surface->isInside(pos));
				sVirtual[i]->translate(pos);
				sVirtual[i]->rotate(angle);
				sOverlapped[i] = this->surface->check(sVirtual[i]);
				if (sOverlapped[i]!=NULL){
					delete sVirtual[i];
					aVoxels[i]->miss();

					#pragma omp atomic
					missCounter++;

					if(this->getFactor()>FACTOR_LIMIT && this->voxels->analyzeVoxel(aVoxels[i], sOverlapped[i], this->surface)){
						#pragma omp critical
						this->voxels->remove(aVoxels[i]);
					}
				}
			}

//			#pragma omp critical
//			std::cout << "[" << omp_get_thread_num() << "] end looping" << std::endl << std::flush;

		}// parallel

		// processing potentially non overlaping shapes
		for(int i=0; i<tmpSplit; i++){

			t += this->getFactor() * dt;

			if (sOverlapped[i]==NULL){
				if(this->surface->check(sVirtual[i])==NULL){ // next checking of overlapping

					l++;
					sVirtual[i]->no = l;
					sVirtual[i]->time = t;


					if (this->params->modifiedRSA){
						this->modifiedRSA(sVirtual[i], aVoxels[i]);
					}


					if (aVoxels[i]!=this->voxels->getVoxel(aVoxels[i]->getPosition(), aVoxels[i]->getOrientation())){
						Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v1 = this->voxels->getVoxel(aVoxels[i]->getPosition(), aVoxels[i]->getOrientation());
						std::cout << "Problem: PackingGenerator - inconsistent voxels positions: " <<
								" (" << aVoxels[i]->getPosition()[0] << ", " << aVoxels[i]->getPosition()[1] << ")" <<
								", (" << v1->getPosition()[0] << ", " << v1->getPosition()[1] << ")" <<
								std::endl;
					}

					this->surface->add(sVirtual[i]);
					this->packing.push_back(sVirtual[i]);

					if(this->getFactor()>FACTOR_LIMIT){
						this->analyzeRegion(aVoxels[i]);
					}else{
						this->voxels->remove(aVoxels[i]);
						this->voxels->removeTopLevelVoxel(aVoxels[i]);
					}
					missCounter = 0;
				}else{
					delete sVirtual[i];
				}
			}
		}

		//whether splitting voxels
		if (missCounter > 0) { // v.getMissCounter() % iSplit == 0){ //
			missCounter = 0;
			int v0 = this->voxels->length();

			std::cout << "[" << this->seed << " Surface::doIteration] splitting " << v0 << " voxels ";
			std::cout.flush();
//						this->toPovray("snapshot_before_" + std::to_string(snapshotCounter++) + ".pov");

			bool b = voxels->splitVoxels(this->params->minDx, this->params->maxVoxels, this->surface->getNeighbourGrid(), this->surface);
			int v1 = this->voxels->length();

//					this->toPovray("snapshot_after_" + std::to_string(snapshotCounter++) + ".pov");
			std::cout << " done: " << this->packing.size() << " shapes, " << v1 << " voxels, new voxel size " << voxels->getVoxelSize() << ", factor " << this->getFactor() << std::endl;

			if (b) {
				tmpSplit *=  ((double)v1 / v0);
				if (tmpSplit > 0.5 * std::numeric_limits<int>::max())
					tmpSplit = 0.5 * std::numeric_limits<int>::max();
				if(voxels->length()<1000 && tmpSplit>100)
					tmpSplit /= 10.0;

				delete[] sOverlapped;
				delete[] sVirtual;
				delete[] aVoxels;
				sOverlapped = new Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *[tmpSplit];
				sVirtual = new Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *[tmpSplit];
				aVoxels = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *[tmpSplit];

			} else {
				this->analyzeVoxels();
			}
		}

//					this->toPovray("test_" + std::to_string(this->voxels->getVoxelSize()) + ".pov");
	} // while

	for(int i=0; i<aRNDSize; i++){
		delete aRND[i];
	}
	delete[] aRND;

	delete this->surface;

	std::cout << "[" << seed << " PackingGenerator::createPacking] finished after generating " << l << " shapes" << std::endl;
}

#else

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::createPacking(){

	std::cout << "[" << this->seed << " PackingGenerator::createPacking] started" << std::endl;
	#ifdef _OPENMP
	std::cout << "[" << this->seed << " PackingGenerator::createPacking] using up to " << omp_get_max_threads() << " concurent treads" << std::endl;
	#endif
	int missCounter = 0;
	RND rnd(this->seed);
	ShapeFactory::initShapeClass(this->params->particleType, this->params->particleAttributes);
	Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s = ShapeFactory::createShape(&rnd);

	this->voxels = new VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>(this->params->surfaceSize, s->getVoxelSize());

	double gridSize = s->getNeighbourListCellSize();
	if (gridSize < this->params->thresholdDistance)
		gridSize = this->params->thresholdDistance;


	if(params->boundaryConditions.compare("free")==0)
		this->surface = new NBoxFBC(SPATIAL_DIMENSION, this->params->surfaceSize, s->getNeighbourListCellSize(), s->getVoxelSize());
	else
		this->surface = new NBoxPBC(SPATIAL_DIMENSION, this->params->surfaceSize, s->getNeighbourListCellSize(), s->getVoxelSize());

	double dt = s->getVolume() / this->surface->getArea();
	delete s;
	int l = 0;
	double t = 0;
	int tmpSplit = this->params->split;
//	int snapshotCounter = 0;


	while (!this->isSaturated() && t<params->maxTime && missCounter<params->maxTriesWithoutSuccess) {
		t += this->getFactor() * dt;
		s = ShapeFactory::createShape(&rnd);

		Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v;
		double pos[SPATIAL_DIMENSION];
		double angle[ANGULAR_DIMENSION];
		do{
			v = this->voxels->getRandomVoxel(&rnd);
			this->voxels->getRandomPositionAndOrientation(pos, angle, v, &rnd);
		}while(!this->surface->isInside(pos));
		s->translate(pos);
		s->rotate(angle);

		Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *sTmp = this->surface->check(s);
		if (sTmp==NULL) { // if no overlap detected
			l++;
			s->no = l;
			s->time = t;


			if (this->params->modifiedRSA){
				this->modifiedRSA(s, v);
			}


			if (v!=this->voxels->getVoxel(v->getPosition(), v->getOrientation())){
				Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v1 = this->voxels->getVoxel(v->getPosition(), v->getOrientation());
				std::cout << "Problem: PackingGenerator - inconsistent voxels positions: " <<
						" (" << v->getPosition()[0] << ", " << v->getPosition()[1] << ")" <<
						", (" << v1->getPosition()[0] << ", " << v1->getPosition()[1] << ")" <<
						std::endl;
			}

			this->surface->add(s);
			this->packing.push_back(s);

			if(this->getFactor()>FACTOR_LIMIT){
				this->analyzeRegion(v);
			}else{
				this->voxels->remove(v);
				this->voxels->removeTopLevelVoxel(v);
			}
			missCounter = 0;
//			std::cout << ((Ellipse *)s)->toWolfram() << std::endl;
		}else{ // overlap detected
			v->miss();
			missCounter++;

			if(this->getFactor()>FACTOR_LIMIT && this->voxels->analyzeVoxel(v, sTmp, this->surface))
				this->voxels->remove(v);
			else if (v->getMissCounter() % this->params->analyze == 0) {
//				if(this->voxels->analyzeVoxel(v, this->surface->getNeighbourGrid(), this->surface) && this->getFactor()>FACTOR_LIMIT)
//					this->analyzeRegion(v);
			}
			if (missCounter > tmpSplit) { // v.getMissCounter() % iSplit == 0){ //
				missCounter = 0;
				int v0 = this->voxels->length();

				std::cout << "[" << this->seed << " Surface::doIteration] splitting " << v0 << " voxels ";
				std::cout.flush();
//				this->toPovray("snapshot_before_" + std::to_string(snapshotCounter++) + ".pov");

				bool b = voxels->splitVoxels(this->params->minDx, this->params->maxVoxels, this->surface->getNeighbourGrid(), this->surface);
				int v1 = this->voxels->length();

//				this->toPovray("snapshot_after_" + std::to_string(snapshotCounter++) + ".pov");
				std::cout << " done: " << this->packing.size() << " shapes, " << v1 << " voxels, new voxel size " << voxels->getVoxelSize() << ", factor " << this->getFactor() << std::endl;

				if (b) {
					tmpSplit *=  ((double)v1 / v0);
					if (tmpSplit > 0.5 * std::numeric_limits<int>::max())
						tmpSplit = 0.5 * std::numeric_limits<int>::max();
					if(voxels->length()<1000 && tmpSplit>100)
						tmpSplit /= 10.0;
				} else {
					this->analyzeVoxels();
				}
//				this->toPovray("test_" + std::to_string(this->voxels->getVoxelSize()) + ".pov");
//				this->toWolfram("test_" + std::to_string(this->voxels->getVoxelSize()) + ".nb");
			}
			delete s;
		} // else
	} // while
	delete this->surface;

	std::cout << "[" << seed << " PackingGenerator::createPacking] finished after generating " << l << " shapes" << std::endl;
}
#endif

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::run(){
	this->createPacking();
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> * PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::getPacking(){
	return &this->packing;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toPovray(std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> * packing, double size, VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *voxels, std::string filename){
	std::ofstream file(filename);

	file << "#include \"colors.inc\"" << std::endl;
	file << "background { color White }" << std::endl;
	file << "camera { orthographic location <" << size / 2 << ", " << size / 2 << ", " << (1.3 * size) << "> look_at  <" << size / 2 << ", " << size / 2 << ",  0> }" << std::endl;
	file << "light_source { < 1000.0, 1000.0, 1000.0> color White shadowless parallel point_at <" << size / 2 << ", " << size / 2 << ",  0>}" << std::endl;
	file << "#declare layer=union{" << std::endl;

	file << "  polygon {5, <0.0, 0.0, 0.0>, <0.0, " << size << ", 0.0>, <" << size << ", " << size << ", 0.0>, <" << size << ", 0.0, 0.0>, <0.0, 0.0, 0.0>  texture { finish { ambient 1 diffuse 0 } pigment { color Gray} } }" << std::endl;
//	file << "  text { ttf \"timrom.ttf\" \"0\" 1, 0 pigment { color Black } scale 1.0 translate < 0, 0, 0.0002> }" << std::endl;

	for (Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s : *packing) {
		double *da = s->getPosition();
		file << "  text { ttf \"timrom.ttf\" \"" << s->no << "\" 1, 0 pigment { color White } scale 0.2 translate < " << da[0] << ", " << da[1] << ", 0.01> }" << std::endl;
		file << s->toPovray();
	}

	file << "}" << std::endl;

	if(voxels!=NULL){
		file << "#declare voxels=union{" << std::endl;
		file << voxels->toPovray() << std::endl;
		file << "}" << std::endl;
	}
	file << "#declare result=union{" << std::endl;
	file << "  object { layer }" << std::endl;

	if (voxels!=NULL)
		file << "  object { voxels }" << std::endl;
	file << "}" << std::endl;
	file << "object{ result	rotate x*360*clock }" << std::endl;

	file.close();
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toPovray(std::string filename){
	PackingGenerator::toPovray(&(this->packing), this->params->surfaceSize, this->voxels, filename);
}


template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toWolfram(std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> * packing, double size, VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *voxels, std::string filename){
	std::ofstream file(filename);

	file << "Graphics[{Red";

//	file << "  polygon {4, <0, 0, 0.0>, <0, " << size << ", 0.0>, <" << size << ", " << size << ", 0.0>, <" << size << ", 0, 0.0>  texture { finish { ambient 1 diffuse 0 } pigment { color Gray} } }" << std::endl;
//	file << "  text { ttf \"timrom.ttf\" \"0\" 1, 0 pigment { color Black } scale 1.0 translate < 0, 0, 0.0002> }" << std::endl;

	for (Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s : *packing) {
//		double *da = s->getPosition();
//		file << "  text { ttf \"timrom.ttf\" \"" << s->no << "\" 1, 0 pigment { color White } scale 0.2 translate < " << da[0] << ", " << da[1] << ", 5> }" << std::endl;
		file << ", " << std::endl << s->toWolfram();
	}


	if(voxels!=NULL){
		file << ", Black, " << std::endl << voxels->toWolfram();
	}

	file << std::endl << "}]" << std::endl;

	file.close();
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toWolfram(std::string filename){
	PackingGenerator::toWolfram(&(this->packing), this->params->surfaceSize, this->voxels, filename);
}

