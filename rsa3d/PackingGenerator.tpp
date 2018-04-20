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
#include <limits>
#include <iomanip>


#ifdef _OPENMP
#include <omp.h>
#endif

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
double PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::FACTOR_LIMIT = 5.0;

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::PackingGenerator(int seed, Parameters *params) {
	this->params = params;
	this->seed = seed;
	RND rnd(this->seed);
	Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s = ShapeFactory::createShape(&rnd);

	this->spatialSize = this->params->surfaceSize;
	this->angularSize = s->getVoxelAngularSize();

	this->voxels = new VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION>(this->spatialSize, s->getVoxelSpatialSize(), this->angularSize);

	double gridSize = s->getNeighbourListCellSize();
	if (gridSize < this->params->thresholdDistance)
		gridSize = this->params->thresholdDistance;


	if(this->params->boundaryConditions.compare("free")==0)
		this->surface = new NBoxFBC(SPATIAL_DIMENSION, this->params->surfaceSize, gridSize, s->getVoxelSpatialSize());
	else
		this->surface = new NBoxPBC(SPATIAL_DIMENSION, this->params->surfaceSize, gridSize, s->getVoxelSpatialSize());

	delete s;
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::~PackingGenerator() {
	for(Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s : this->packing)
		delete s;
	if (this->voxels!=NULL)
		delete this->voxels;
	delete this->surface;
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
int PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::analyzeVoxels(unsigned short depth) {

	std::cout << "[" << this->seed << " PackingGenerator::analyzeVoxels] " << this->voxels->length() << " voxels, " << this->packing.size() << " shapes, factor = " << this->getFactor() << ", depth = " << depth << " " << std::flush;

	int begin = this->voxels->length();
	std::vector< Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION>* > toRemove;

	#pragma omp parallel for
	for (int i = 0; i < this->voxels->length(); i++) {
		Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v = this->voxels->get(i);
		bool bRemove = this->voxels->analyzeVoxel(v, this->surface->getNeighbourGrid(), this->surface, depth);
		if (bRemove) {
			#pragma omp critical
			toRemove.push_back(v);
		}
		if (i%10000 == 0){ std::cout << "."; std::cout.flush(); }
	}
	std::cout << " removing " << toRemove.size() << " voxels "<< std::flush;

	this->voxels->remove(&toRemove);

	std::cout << " done: " << this->voxels->length() << " voxels remained, factor = " << this->getFactor() << std::endl << std::flush;

	return begin - this->voxels->length();
}

#else

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
int PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::analyzeVoxels(unsigned short depth) {

	std::cout << "[" << this->seed << " PackingGenerator::analyzeVoxels] " << this->voxels->length() << " voxels, " << this->packing.size() << " shapes, factor = " << this->getFactor() << " ";

	int begin = this->voxels->length();

	for (int i = 0; i < this->voxels->length(); i++) {
		Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v = this->voxels->get(i);
		bool bRemove = this->voxels->analyzeVoxel(v, this->surface->getNeighbourGrid(), this->surface, depth);
		if (bRemove) {
			this->voxels->remove(v);
			i--;
		}
		if (i%10000 == 0){ std::cout << "."; std::cout.flush(); }
	}

	std::cout << " done: " << this->voxels->length() << " voxels remained, factor = " << this->getFactor() << std::endl;

	return begin - this->voxels->length();
}

#endif

// analyzes all voxels inside a region around v
template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
int PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::analyzeRegion(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v){
	int begin = this->voxels->length();
	std::vector<Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> region;
	this->voxels->getNeighbours(&region, v->getPosition());
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

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
bool PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::isInside(double *position, double *orientation){
	for(unsigned short i=0; i<SPATIAL_DIMENSION; i++){
		if (position[i]>=this->spatialSize || position[i]<0)
			return false;
	}
	for(unsigned short i=0; i<ANGULAR_DIMENSION; i++){
		if (orientation[i]>=this->angularSize || orientation[i]<0)
			return false;
	}
	return true;

}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::testPacking(std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> *vShapes, double maxTime){

#ifdef _OPENMP
	int maxthreads = omp_get_max_threads();
#else
	int maxthreads = 1;
#endif

	int loop = 1000*maxthreads;

	std::cout.precision(std::numeric_limits< double >::max_digits10);

	std::cout << "[" << this->seed << " PackingGenerator::testPacking] using up to " << maxthreads << " concurrent treads" << std::endl;


	RND rnd(this->seed);
	Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s = ShapeFactory::createShape(&rnd);
	double dt = s->getVolume() / this->surface->getArea();
	delete s;

	for(Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s : *vShapes){
		this->surface->add(s);
		this->packing.push_back(s);
	}
	std::cout << "[" << this->seed << " PackingGenerator::testPacking] " << vShapes->size() << " shapes restored" << std::endl;



	double t = 0;

	RND **aRND = new RND*[maxthreads];
	for(int i=0; i<maxthreads; i++){
		aRND[i] = new RND((int)(1000*rnd.nextValue()));
	}

	while (t<maxTime) {

		std::cout << "\r" << "[" << this->seed << " PackingGenerator::testPacking] t=" << std::setprecision(4) << t/maxTime << " choosing " << loop << " shapes..." << std::flush;

		#pragma omp parallel for
		for(int i = 0; i<loop; i++){

		#ifdef _OPENMP
			int tid = omp_get_thread_num();
		#else
			int tid = 0;
		#endif
			Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *sVirtual = ShapeFactory::createShape(aRND[tid]);
			Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v;
			double pos[SPATIAL_DIMENSION];
			std::array <double, ANGULAR_DIMENSION> angle;
			do{
				v = this->voxels->getRandomVoxel(aRND[tid]);
				this->voxels->getRandomPositionAndOrientation(pos, angle.data(), v, aRND[tid]);
			}while(!this->isInside(pos, angle.data()));
			// setting shape position and orientation
			sVirtual->translate(pos);
			sVirtual->rotate(angle.data());
			// checking if shape overlaps with any shape in the packing
			if (this->surface->check(sVirtual)==NULL){
				#pragma omp critical(stdout)
				{
					std::cout << std::endl << "\t non overlapping shape found " << std::setprecision(10) << sVirtual->toString() << std::endl << std::flush;
					std::cout << "\t povray: " << std::endl << std::setprecision(10) << sVirtual->toPovray() << std::endl << std::flush;
				}
				double *position = sVirtual->getPosition();
				const double *tmp = sVirtual->getOrientation();
				double orientation[ANGULAR_DIMENSION];
				std::copy(tmp, tmp + ANGULAR_DIMENSION, orientation);
				double delta = 0.0001;
				for(unsigned short j = 0; j< ANGULAR_DIMENSION; j++)
					orientation[j] -= 0.5*delta;

				Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *sCovers = NULL;
				std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> vNeighbours;
				this->surface->getNeighbours(&vNeighbours, position);
				for(Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *sTmp : vNeighbours){
					if (sTmp->pointInside(this->surface, position, orientation, delta)){
						sCovers = sTmp;
						break;
					}
				}
				if (sCovers!=NULL)
				#pragma omp critical(stdout)
				std::cout << "\t in exclusion zone of " << sCovers->toString() << std::endl;
			}
			delete sVirtual;
		} // parallel for

		t += dt * loop;
	} // while

	for(int i=0; i<maxthreads; i++){
		delete aRND[i];
	}
	delete[] aRND;

	std::cout << "[" << seed << " PackingGenerator::testPacking] finished after time " << t << std::endl;
}




#ifdef _OPENMP

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::createPacking(){

	int maxthreads = omp_get_max_threads();

	std::cout.precision(std::numeric_limits< double >::max_digits10);

	std::cout << "[" << this->seed << " PackingGenerator::createPacking] using up to " << omp_get_max_threads() << " concurrent treads" << std::endl;

	int checkedAgain = 0;
	int added = 0;
	int missCounter = 0;
	unsigned short depthAnalyze = 0;


	RND rnd(this->seed);
	Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s = ShapeFactory::createShape(&rnd);
	double dt = s->getVolume() / this->surface->getArea();
	delete s;

	int l = 0;
	double t = 0, factor = 1;
	int tmpSplit = this->params->split, oldTmpSplit = tmpSplit;
//	int snapshotCounter = 0;
	RND **aRND = new RND*[maxthreads];
	for(int i=0; i<maxthreads; i++){
		aRND[i] = new RND((int)(1000*rnd.nextValue()));
	}

	Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> **sOverlapped = new Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *[tmpSplit];
	Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> **sVirtual = new Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *[tmpSplit];
	Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> **aVoxels = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *[tmpSplit];

	while (!this->isSaturated() && t<params->maxTime && missCounter<params->maxTriesWithoutSuccess) {

		std::cout << "[" << this->seed << " PackingGenerator::createPacking] choosing " << tmpSplit << " shapes..." << std::flush;
		factor = this->getFactor();

		#pragma omp parallel for
		for(int i = 0; i<tmpSplit; i++){
			int tid = omp_get_thread_num();

			sVirtual[i] = ShapeFactory::createShape(aRND[tid]);
			double pos[SPATIAL_DIMENSION];
			std::array <double, ANGULAR_DIMENSION> angle;
			do{
				aVoxels[i] = this->voxels->getRandomVoxel(aRND[tid]);
				this->voxels->getRandomPositionAndOrientation(pos, angle.data(), aVoxels[i], aRND[tid]);
			}while(!this->isInside(pos, angle.data()));
			// setting shape position and orientation
			sVirtual[i]->translate(pos);
			sVirtual[i]->rotate(angle.data());
			// checking if shape overlaps with any shape in the packing
			sOverlapped[i] = this->surface->check(sVirtual[i]);
		} // parallel for


		std::cout << " processing shapes..." << std::flush;

		checkedAgain = 0;
		added = 0;

/*
		if(this->voxels->length()<20){
			std::string filename = "shapes_" + std::to_string(this->voxels->length()) + ".dbg";
			std::ofstream file(filename, std::ios::binary);
			file.write((char *)(&tmpSplit), sizeof(int));
			for(int i=0; i<tmpSplit; i++){
				sVirtual[i]->store(file);
			}
			this->store(file);
			file.close();
		}
*/

		// sequentially processing potentially non overlaping shapes
		for(int i=0; i<tmpSplit; i++){

			t += factor * dt;

			// if there were no intersecting particles in the packing
			if (sOverlapped[i]==NULL){
				checkedAgain++;
				// checking overlapping again
				sOverlapped[i] = this->surface->check(sVirtual[i]);

				// if there is still no overlap sVirtula is added to the packing and corresponding voxel is removed
				if(sOverlapped[i]==NULL){ // if non overlapping

					added++;
					l++;
					sVirtual[i]->no = l;
					sVirtual[i]->time = t;


					if (this->params->modifiedRSA){
						this->modifiedRSA(sVirtual[i], aVoxels[i]);
					}

					// consistency check
					if (aVoxels[i]!=this->voxels->getVoxel(aVoxels[i]->getPosition(), aVoxels[i]->getOrientation())){
						Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v = this->voxels->getVoxel(aVoxels[i]->getPosition(), aVoxels[i]->getOrientation());
						std::cout << std::endl << "Problem: PackingGenerator - inconsistent voxels positions: [" << aVoxels[i]->toString() << "], [" << v->toString() << "]" << std::endl;
						std::cout << "size: " << this->voxels->getVoxelSize() << ", angular size: " << this->voxels->getVoxelAngularSize() << std::endl;
						std::cout << "shape: " << sVirtual[i]->toString() << std::endl;

					}

					this->surface->add(sVirtual[i]);
					this->packing.push_back(sVirtual[i]);
					this->voxels->remove(aVoxels[i]);
					this->voxels->removeTopLevelVoxel(aVoxels[i]);
				}  //second overlapping check
			} // first overlapping check
		}  // for


		// now processing voxels, with missed hits

		std::cout << " processing voxels..." << std::flush;


		// when factor is large there is no need to check voxels because the chance of finding removing voxel is very small
		if (factor<1000){
			#pragma omp parallel for
			for(int i=0; i<tmpSplit; i++){
				// non overlapping shapes was already processed. Now we process missed hits (where there was an overlap)
				if (sOverlapped[i]!=NULL){
					Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *vTmp = this->voxels->getVoxel(sVirtual[i]->getPosition(), sVirtual[i]->getOrientation());
					if (vTmp != NULL){
						// checking voxels consistency
						if (vTmp != aVoxels[i]){
							std::cout << std::endl << "Problem: PackingGenerator - inconsistent voxels: " << i << std::flush;
						}
						vTmp->miss();
						if( (vTmp->getMissCounter() > 0 &&  vTmp->getMissCounter() % this->params->analyze == 0) ){
							if (this->voxels->analyzeVoxel(vTmp, sOverlapped[i], this->surface)){
								#pragma omp critical
								this->voxels->remove(vTmp);
							}
						}
					}
					delete sVirtual[i];
				}
			} // parallel for
		}

		// processing remaining voxels with miss

		std::cout << "done, double checked: " << checkedAgain << " added: " << added << ", time: " << t << ", shapes: " << l << std::endl << std::flush;

		//whether splitting voxels
		if (added == 0) { // v.getMissCounter() % iSplit == 0){ //
			missCounter += tmpSplit;
			int v0 = this->voxels->length(), v1 = v0;

			std::cout << "[" << this->seed << " PackingGenerator::createPacking] splitting " << v0 << " voxels ";
			std::cout.flush();
//						this->toPovray("snapshot_before_" + std::to_string(snapshotCounter++) + ".pov");

			bool b = voxels->splitVoxels(this->params->minDx, this->params->maxVoxels, this->surface->getNeighbourGrid(), this->surface);
			if (b){
				v1 = this->voxels->length();
//				this->toPovray("snapshot_after_" + std::to_string(snapshotCounter++) + ".pov");
				std::cout << " done. " << this->packing.size() << " shapes, " << v1 << " voxels, new voxel size: " << voxels->getVoxelSize() << ", angular size: " << this->voxels->getVoxelAngularSize() << ", factor: " << this->getFactor() << std::endl;
				missCounter = 0;
			}else{
				std::cout << "skipped." << std::endl;
				this->analyzeVoxels(depthAnalyze);
				tmpSplit = 1.1*tmpSplit;
				v1 = this->voxels->length();
			}
			// if number of voxels has changed
			if (v1!=v0){
				tmpSplit *= ((double)v1 / v0);
			}else{
				tmpSplit = 1.1*tmpSplit + omp_get_max_threads();
			}

			if (tmpSplit > std::max(this->params->maxVoxels/20, 10*this->params->split))
				tmpSplit = std::max(this->params->maxVoxels/20, 10*this->params->split);
			if(voxels->length()<0.001*this->params->maxVoxels && tmpSplit > 10*omp_get_max_threads())
				tmpSplit /= 10.0;
			if(tmpSplit < 10*omp_get_max_threads())
				tmpSplit = 10*omp_get_max_threads();

			if (!b && (double)(v0-v1)/v0 < 0.1){ // not much voxels removed
				depthAnalyze++;
			}else{
				if (depthAnalyze>0)
					depthAnalyze--;
			}

			if(tmpSplit != oldTmpSplit){

				delete[] sOverlapped;
				delete[] sVirtual;
				delete[] aVoxels;

				sOverlapped = new Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *[tmpSplit];
				sVirtual = new Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *[tmpSplit];
				aVoxels = new Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *[tmpSplit];

				oldTmpSplit = tmpSplit;
			}
//			this->printRemainingVoxels("voxels_" + std::to_string(this->voxels->getVoxelSize()));
//			this->toWolfram("test_" + std::to_string(this->voxels->getVoxelSize()) + ".nb");
//			this->toPovray("test_" + std::to_string(this->voxels->getVoxelSize()) + ".pov");
//			std::string filename = "snapshot_" + std:::to_string(this->packing.size()) + "_" + std::to_string(this->voxels->length()) + ".dbg";
//			std::ofstream file(filename, std::ios::binary);
//			this->store(file);
//			file.close();
		}else{
			missCounter = 0;
			if (depthAnalyze>0)
				depthAnalyze--;
		}
	} // while

	for(int i=0; i<maxthreads; i++){
		delete aRND[i];
	}
	delete[] aRND;

	std::cout << "[" << seed << " PackingGenerator::createPacking] finished after generating " << l << " shapes" << std::endl;
}

#else

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::createPacking(){
	std::cout.precision(std::numeric_limits< double >::max_digits10);

	std::cout << "[" << this->seed << " PackingGenerator::createPacking] started" << std::endl;
	int missCounter = 0;

	RND rnd(this->seed);
	Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s = ShapeFactory::createShape(&rnd);
	double dt = s->getVolume() / this->surface->getArea();
	delete s;

	int l = 0;
	double t = 0;
	int tmpSplit = this->params->split;
	int depthAnalyze = 0;
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
		}while(!this->isInside(pos, angle));
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
				std::cout << std::endl << "Problem: PackingGenerator - inconsistent voxels positions: [" << v->toString()<< "], [" << v1->toString() << "]" << std::endl;
				std::cout << "size: " << this->voxels->getVoxelSize() << ", angular size: " << this->voxels->getVoxelAngularSize() << std::endl;
				std::cout << "shape: " << s->toString() << std::endl;
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
			if (depthAnalyze>0)
				depthAnalyze--;
//			std::cout << ((Ellipse *)s)->toWolfram() << std::endl;
		}else{ // overlap detected
			v->miss();
			missCounter++;

			if(this->getFactor()>FACTOR_LIMIT && (v->getMissCounter() > 0 &&  v->getMissCounter() % this->params->analyze == 0) && this->voxels->analyzeVoxel(v, sTmp, this->surface)){
				this->voxels->remove(v);
			}
			if (missCounter > tmpSplit) { // v.getMissCounter() % iSplit == 0){ //
				missCounter = 0;
				int v0 = this->voxels->length(), v1 = v0;

				std::cout << "[" << this->seed << " PackingGenerator::createPacking] splitting " << v0 << " voxels ";
				std::cout.flush();
//				this->toPovray("snapshot_before_" + std::to_string(snapshotCounter++) + ".pov");

				bool b = voxels->splitVoxels(this->params->minDx, this->params->maxVoxels, this->surface->getNeighbourGrid(), this->surface);
				if (b){
					v1 = this->voxels->length();
	//				this->toPovray("snapshot_after_" + std::to_string(snapshotCounter++) + ".pov");
					std::cout << " done. " << this->packing.size() << " shapes, " << v1 << " voxels, new voxel size: " << voxels->getVoxelSize() << ", angular size: " << this->voxels->getVoxelAngularSize() << ", factor: " << this->getFactor() << std::endl;
					missCounter = 0;
				}else{
					std::cout << "skipped." << std::endl;
					this->analyzeVoxels(depthAnalyze);
					tmpSplit = 1.1*tmpSplit;
					v1 = this->voxels->length();
				}
				// if number of voxels has changed
				if (v1!=v0){
					tmpSplit *= ((double)v1 / v0);
				}else{
					tmpSplit = 1.1*tmpSplit + 10;
				}

				if (tmpSplit > std::max(this->params->maxVoxels/20, 10*this->params->split))
					tmpSplit = std::max(this->params->maxVoxels/20, 10*this->params->split);
				if(voxels->length()<0.001*this->params->maxVoxels && tmpSplit > 100)
					tmpSplit /= 10.0;
				if(tmpSplit < 100)
					tmpSplit = 100;

				if (!b && v0==v1){ // if nothing changed with voxels
					depthAnalyze++;
				}else{
					if (depthAnalyze>0)
						depthAnalyze--;
				}
//				this->printRemainingVoxels("voxels_" + std::to_string(this->voxels->getVoxelSize()));
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
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toPovray(std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> * packing, double size, VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *voxels, const std::string &filename){
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
// periodic boundary check
		double tr[2];

		if (da[0]<0.1*size){
			tr[0] = size;
			tr[1] = 0;
			s->translate(tr);
			file << s->toPovray();
			tr[0] = -tr[0];
			tr[1] = -tr[1];
			s->translate(tr);
			if (da[1]<0.1*size){
				tr[0] = size;
				tr[1] = size;
				s->translate(tr);
				file << s->toPovray();
				tr[0] = -tr[0];
				tr[1] = -tr[1];
				s->translate(tr);
			}
			if(da[1]>0.9*size){
				tr[0] = size;
				tr[1] = -size;
				s->translate(tr);
				file << s->toPovray();
				tr[0] = -tr[0];
				tr[1] = -tr[1];
				s->translate(tr);
			}
		}
		if (da[0]>0.9*size){
			tr[0] = -size;
			tr[1] = 0;
			s->translate(tr);
			file << s->toPovray();
			tr[0] = -tr[0];
			tr[1] = -tr[1];
			s->translate(tr);
			if (da[1]<0.1*size){
				tr[0] = -size;
				tr[1] = size;
				s->translate(tr);
				file << s->toPovray();
				tr[0] = -tr[0];
				tr[1] = -tr[1];
				s->translate(tr);
			}
			if(da[1]>0.9*size){
				tr[0] = -size;
				tr[1] = -size;
				s->translate(tr);
				file << s->toPovray();
				tr[0] = -tr[0];
				tr[1] = -tr[1];
				s->translate(tr);
			}
		}
		if (da[1]<0.1*size){
			tr[0] = 0;
			tr[1] = size;
			s->translate(tr);
			file << s->toPovray();
			tr[0] = -tr[0];
			tr[1] = -tr[1];
			s->translate(tr);
		}
		if (da[1]>0.9*size){
			tr[0] = 0;
			tr[1] = -size;
			s->translate(tr);
			file << s->toPovray();
			tr[0] = -tr[0];
			tr[1] = -tr[1];
			s->translate(tr);
		}

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
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toPovray(const std::string &filename){
	PackingGenerator::toPovray(&(this->packing), this->params->surfaceSize, this->voxels, filename);
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toWolfram(double *da, const std::string &filename){

	std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> vShapes;
	this->surface->getNeighbours(&vShapes, da);

	std::vector<Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> vVoxels;
	this->voxels->getNeighbours(&vVoxels, da);

	std::ofstream file(filename);
	file << "Graphics[{Red";

	for (Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s : vShapes) {
		file << ", " << std::endl << s->toWolfram();
	}
	if (vVoxels.size()>0){
		file << ", Black, " << std::endl;
		for(Voxel<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *v: vVoxels){
			file << v->toWolfram(this->voxels->getVoxelSize(), this->voxels->getVoxelAngularSize());
		}
	}
	file << std::endl << "}]" << std::endl;
	file.close();
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::printRemainingVoxels(const std::string &prefix){
	if (this->voxels->length()>20)
		return;
	for(int i=0; i<this->voxels->length(); i++){
		std::string filename(prefix + "_" + std::to_string(i) + ".nb");
		this->toWolfram(this->voxels->getVoxel(i)->getPosition(), filename);
	}
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toWolfram(std::vector<Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *> * packing, double size, VoxelList<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *voxels, const std::string &filename){
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
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toWolfram(const std::string &filename){
	PackingGenerator::toWolfram(&(this->packing), this->params->surfaceSize, this->voxels, filename);
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::toFile(const std::string &filename) {
	std::ofstream file(filename, std::ios::binary);
    if (!file)
        die("Cannot open file " + filename + " to store packing");
	for (Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *s : this->packing) {
		s->store(file);
	}
	file.close();
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::store(std::ostream &f) const{
	unsigned short sd = SPATIAL_DIMENSION;
	unsigned short ad = ANGULAR_DIMENSION;
	f.write((char *)(&sd), sizeof(unsigned char));
	if (ad>0)
		f.write((char *)(&ad), sizeof(unsigned char));
	int size = this->packing.size();
	f.write((char *)(&size), sizeof(int));

	for(Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s: this->packing){
		s->store(f);
	}
	this->voxels->store(f);
}

template <unsigned short SPATIAL_DIMENSION, unsigned short ANGULAR_DIMENSION>
void PackingGenerator<SPATIAL_DIMENSION, ANGULAR_DIMENSION>::restore(std::istream &f){
	unsigned char sd = SPATIAL_DIMENSION;
	unsigned char ad = ANGULAR_DIMENSION;

	f.read((char *)(&sd), sizeof(unsigned char));
	if (ad > 0)
		f.read((char *)(&ad), sizeof(unsigned char));

	if (sd!=SPATIAL_DIMENSION || ad!=ANGULAR_DIMENSION){
		std::cout << "[ERROR] cannot restore PackingGenerator: incompatible dimensions: read " << f.gcount() << " bytes." << std::endl;
		return;
	}
	int size;
	RND rnd;
	f.read((char *)(&size), sizeof(int));
	this->packing.clear();
	this->surface->clear();
	for(int i=0; i<size; i++){
		Shape<SPATIAL_DIMENSION, ANGULAR_DIMENSION> *s = ShapeFactory::createShape(&rnd);
		s->restore(f);
		this->surface->add(s);
		this->packing.push_back(s);
//		std::cout << s->toString() << std::endl;
	}
	this->voxels->restore(f);
}
