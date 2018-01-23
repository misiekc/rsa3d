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

template <ushort DIMENSION>
double PackingGenerator<DIMENSION>::FACTOR_LIMIT = 5.0;

template <ushort DIMENSION>
PackingGenerator<DIMENSION>::PackingGenerator(int s, Parameters *p) {
	this->params = p;
	this->seed = s;
	this->voxels = NULL;
	this->surface = NULL;

}

template <ushort DIMENSION>
PackingGenerator<DIMENSION>::~PackingGenerator() {
	for(Shape<DIMENSION> *s : this->packing)
		delete s;
	if (this->voxels!=NULL)
		delete this->voxels;
}

template <ushort DIMENSION>
bool PackingGenerator<DIMENSION>::isSaturated() {
	return (this->voxels->length() == 0);
}

template <ushort DIMENSION>
double PackingGenerator<DIMENSION>::getFactor() {
	return this->surface->getArea() / this->voxels->getVoxelsSurface();
}

template <ushort DIMENSION>
int PackingGenerator<DIMENSION>::analyzeVoxels() {

	std::cout << "[" << this->seed << " PackingGenerator::analyzeVoxels] " << this->voxels->length() << " voxels, " << this->packing.size() << " shapes, factor = " << this->getFactor() << std::endl;

	int begin = this->voxels->length();
	int timestamp = this->packing.size();
	for (int i = 0; i < this->voxels->length(); i++) {
		Voxel<DIMENSION> *v = this->voxels->get(i);
		if (this->voxels->analyzeVoxel(v, this->surface->getNeighbourGrid(), this->surface, timestamp)) {
			this->voxels->remove(v);
			i--;
		}
	}

	std::cout << "[" << this->seed << " Surface::analyzeVoxels] " << this->voxels->length() << " voxels remained, factor = " << this->getFactor() << std::endl;

	return begin - this->voxels->length();
}

// analyzes all voxels inside a region around v
template <ushort DIMENSION>
int PackingGenerator<DIMENSION>::analyzeRegion(Voxel<DIMENSION> *v){
	int begin = this->voxels->length();
	std::vector<Voxel<DIMENSION> *> region;
	this->voxels->getNeighbours(&region, v);
	for(Voxel<DIMENSION> *v1: region){
		if (this->voxels->analyzeVoxel(v1, this->surface->getNeighbourGrid(), this->surface))
			this->voxels->remove(v1);
	}
	return begin - this->voxels->length();
}

template <ushort DIMENSION>
void PackingGenerator<DIMENSION>::modifiedRSA(Shape<DIMENSION> *s, Voxel<DIMENSION> *v){
	double da[DIMENSION];

	Shape<DIMENSION> *sn = (Shape<DIMENSION> *)this->surface->getClosestNeighbour(s->getPosition(), NULL);
	if (sn==NULL)
		sn = (Shape<DIMENSION> *)this->surface->getClosestNeighbour(s->getPosition(), &(this->packing));
	if (sn!=NULL){
		double *spos = s->getPosition();
		double *snpos = sn->getPosition();

		double d = sqrt(this->surface->distance2(spos, snpos));
		if (d < this->params->thresholdDistance){

			for(ushort i=0; i<DIMENSION; i++){
				da[i] = (snpos[i] - spos[i]);
			}
			this->surface->vector(da);
			double mindist = s->minDistance(sn);
			for(ushort i=0; i<DIMENSION; i++){
				da[i] = da[i]/d;
				spos[i] += (d-mindist)*da[i];
			}
			this->surface->checkPosition(spos);
			v = this->voxels->getVoxel(spos);
			if (v==NULL){
				std::cout << "Problem: PackingGenerator - voxel not found: " <<
				" (" << spos[0] << ", " << spos[1] << ")" << std::endl;
				exit(0);
			}
		}
	}
}




//#ifdef _OPENMP
template <ushort DIMENSION>
void PackingGenerator<DIMENSION>::createPacking(){

	std::cout << "[" << this->seed << " PackingGenerator::createPacking] started" << std::endl;

	int missCounter = 0;
	RND rnd(this->seed);
	ShapeFactory::initShapeClass(this->params->particleType, this->params->particleAttributes);
	Shape<DIMENSION> *s = ShapeFactory::createShape(&rnd);

	this->voxels = new VoxelList<DIMENSION>(this->params->surfaceSize, s->getVoxelSize());

	double gridSize = s->getNeighbourListCellSize();
	if (gridSize < this->params->thresholdDistance)
		gridSize = this->params->thresholdDistance;


	if(params->boundaryConditions.compare("free")==0)
		this->surface = new NBoxFBC(DIMENSION, this->params->surfaceSize, s->getNeighbourListCellSize(), s->getVoxelSize());
	else
		this->surface = new NBoxPBC(DIMENSION, this->params->surfaceSize, s->getNeighbourListCellSize(), s->getVoxelSize());

	double dt = s->getVolume() / this->surface->getArea();
	delete s;
	int l = 0;
	double t = 0;
	int tmpSplit = this->params->split;
//	int snapshotCounter = 0;


	while (!this->isSaturated() && t<params->maxTime && missCounter<params->maxTriesWithoutSuccess) {


		Shape<DIMENSION> **sOverlapped = new Shape<DIMENSION> *[tmpSplit];
		Shape<DIMENSION> **sVirtual = new Shape<DIMENSION> *[tmpSplit];
		Voxel<DIMENSION> **aVoxels = new Voxel<DIMENSION> *[tmpSplit];

		for(int i = 0; i<tmpSplit; i++){
			sVirtual[i] = ShapeFactory::createShape(&rnd);

			double da[DIMENSION];
			do{
				aVoxels[i] = this->voxels->getRandomVoxel(&rnd);
				this->voxels->getRandomPosition(da, aVoxels[i], &rnd);
			}while(!this->surface->isInside(da));
			sVirtual[i]->translate(da);
		}

		#pragma omp parallel for
		for(int i = 0; i<tmpSplit; i++){
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


					if (aVoxels[i]!=this->voxels->getVoxel(aVoxels[i]->getPosition())){
						Voxel<DIMENSION> *v1 = this->voxels->getVoxel(aVoxels[i]->getPosition());
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

		delete[] sOverlapped;
		delete[] sVirtual;
		delete[] aVoxels;

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
			} else {
				this->analyzeVoxels();
			}
		}
//					this->toPovray("test_" + std::to_string(this->voxels->getVoxelSize()) + ".pov");
	} // while
	delete this->surface;

	std::cout << "[" << seed << " PackingGenerator::createPacking] finished after generating " << l << " shapes" << std::endl;
}

/*
#else

template <ushort DIMENSION>
void PackingGenerator<DIMENSION>::createPacking(){

	std::cout << "[" << this->seed << " PackingGenerator::createPacking] started" << std::endl;

	int missCounter = 0;
	RND rnd(this->seed);
	ShapeFactory::initShapeClass(this->params->particleType, this->params->particleAttributes);
	Shape<DIMENSION> *s = ShapeFactory::createShape(&rnd);

	this->voxels = new VoxelList<DIMENSION>(this->params->surfaceSize, s->getVoxelSize());

	double gridSize = s->getNeighbourListCellSize();
	if (gridSize < this->params->thresholdDistance)
		gridSize = this->params->thresholdDistance;


	if(params->boundaryConditions.compare("free")==0)
		this->surface = new NBoxFBC(DIMENSION, this->params->surfaceSize, s->getNeighbourListCellSize(), s->getVoxelSize());
	else
		this->surface = new NBoxPBC(DIMENSION, this->params->surfaceSize, s->getNeighbourListCellSize(), s->getVoxelSize());

	double dt = s->getVolume() / this->surface->getArea();
	delete s;
	int l = 0;
	double t = 0;
	int tmpSplit = this->params->split;
//	int snapshotCounter = 0;


	while (!this->isSaturated() && t<params->maxTime && missCounter<params->maxTriesWithoutSuccess) {
		t += this->getFactor() * dt;
		s = ShapeFactory::createShape(&rnd);

		Voxel<DIMENSION> *v;
		double da[DIMENSION];
		do{
			v = this->voxels->getRandomVoxel(&rnd);
			this->voxels->getRandomPosition(da, v, &rnd);
		}while(!this->surface->isInside(da));
		s->translate(this->voxels->getRandomPosition(da, v, &rnd));

		Shape<DIMENSION> *sTmp = this->surface->check(s);
		if (sTmp==NULL) { // if no overlap detected
			l++;
			s->no = l;
			s->time = t;


			if (this->params->modifiedRSA){
				this->modifiedRSA(s, v);
			}


			if (v!=this->voxels->getVoxel(v->getPosition())){
				Voxel<DIMENSION> *v1 = this->voxels->getVoxel(v->getPosition());
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
			}
			delete s;
		} // else
	} // while
	delete this->surface;

	std::cout << "[" << seed << " PackingGenerator::createPacking] finished after generating " << l << " shapes" << std::endl;
}
#endif
*/

template <ushort DIMENSION>
void PackingGenerator<DIMENSION>::run(){
	this->createPacking();
}

template <ushort DIMENSION>
std::vector<Shape<DIMENSION> *> * PackingGenerator<DIMENSION>::getPacking(){
	return &this->packing;
}

template <ushort DIMENSION>
void PackingGenerator<DIMENSION>::toPovray(std::vector<Shape<DIMENSION> *> * packing, double size, VoxelList<DIMENSION> *voxels, std::string filename){
	std::ofstream file(filename);

	file << "#include \"colors.inc\"" << std::endl;
	file << "background { color White }" << std::endl;
	file << "camera { orthographic location <" << size / 2 << ", " << size / 2 << ", " << (1.3 * size) << "> look_at  <" << size / 2 << ", " << size / 2 << ",  0> }" << std::endl;
	file << "light_source { < 1000.0, 1000.0, 1000.0> color White shadowless parallel point_at <" << size / 2 << ", " << size / 2 << ",  0>}" << std::endl;
	file << "#declare layer=union{" << std::endl;

	file << "  polygon {4, <0, 0, 0.0>, <0, " << size << ", 0.0>, <" << size << ", " << size << ", 0.0>, <" << size << ", 0, 0.0>  texture { finish { ambient 1 diffuse 0 } pigment { color Gray} } }" << std::endl;
	file << "  text { ttf \"timrom.ttf\" \"0\" 1, 0 pigment { color Black } scale 1.0 translate < 0, 0, 0.0002> }" << std::endl;

	for (Shape<DIMENSION> *s : *packing) {
//		double *da = s->getPosition();
//		file << "  text { ttf \"timrom.ttf\" \"" << s->no << "\" 1, 0 pigment { color White } scale 0.2 translate < " << da[0] << ", " << da[1] << ", 5> }" << std::endl;
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

template <ushort DIMENSION>
void PackingGenerator<DIMENSION>::toPovray(std::string filename){
	PackingGenerator::toPovray(&(this->packing), this->params->surfaceSize, this->voxels, filename);
}

