/*
 * PackingGenerator.cpp
 *
 *  Created on: 16.04.2017
 *      Author: ciesla
 */

#include <memory>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <dirent.h>

#include "PackingGenerator.h"
#include "surfaces/NBoxPBC.h"
#include "surfaces/NBoxFBC.h"
#include "shape/ShapeFactory.h"
#include "shape/ConvexShape.h"
#include "ThreadLocalRND.h"
#include "utils/OMPMacros.h"


double PackingGenerator::FACTOR_LIMIT = 5.0;


PackingGenerator::PackingGenerator(int seed, const Parameters *params) {
	this->params = *params;
	this->seed = seed;
	RND rnd(this->seed);

	this->spatialSize = this->params.surfaceSize;
	this->angularSize = RSAShape::getAngularVoxelSize();
	if (this->params.requestedAngularVoxelSize > this->angularSize)
		this->params.requestedAngularVoxelSize = this->angularSize;

	this->voxels = ShapeFactory::createVoxelList(params->particleType, this->params.surfaceDimension, this->spatialSize, RSAShape::getVoxelSpatialSize(), this->angularSize, this->params.requestedAngularVoxelSize);

	double gridSize = RSAShape::getNeighbourListCellSize();
	if (gridSize < this->params.thresholdDistance)
		gridSize = this->params.thresholdDistance;


	if (this->params.boundaryConditions == "free")
		this->surface = new NBoxFBC(this->params.surfaceDimension, this->params.surfaceSize, gridSize, RSAShape::getVoxelSpatialSize());
	else
		this->surface = new NBoxPBC(this->params.surfaceDimension, this->params.surfaceSize, gridSize, RSAShape::getVoxelSpatialSize());
}


PackingGenerator::~PackingGenerator() {
	delete this->voxels;
	delete this->surface;
}


bool PackingGenerator::isSaturated() {
	return (this->voxels->getLength() == 0);
}


double PackingGenerator::getFactor() {
	return this->surface->getArea() / this->voxels->getVoxelsVolume();
}

void PackingGenerator::modifiedRSA(RSAShape *s, Voxel *v){

	const RSAShape *sn = this->surface->getClosestNeighbour(s->getPosition());
	if (sn == nullptr)
		sn = this->surface->getClosestNeighbour(s->getPosition(), this->packing.getVector());
	if (sn != nullptr){
		RSAVector spos = s->getPosition();
		RSAVector snpos = sn->getPosition();

		double d = sqrt(this->surface->distance2(spos, snpos));
		if (d < this->params.thresholdDistance){
			RSAVector da = this->surface->vector(snpos - spos) / d;
			double mindist = s->minDistance(sn);
			spos += (d - mindist) * da;
			spos = this->surface->checkPosition(spos);
			v = this->voxels->getVoxel(spos, s->getOrientation());
			if (v == nullptr){
				std::cout << "Problem: PackingGenerator - voxel not found: " <<
				" (" << spos[0] << ", " << spos[1] << ")" << std::endl;
				exit(0);
			}
		}
	}
}


bool PackingGenerator::isInside(const RSAVector &position, RSAOrientation &orientation){
	for(unsigned short i=0; i<RSA_SPATIAL_DIMENSION; i++){
		if (position[i]>=this->spatialSize || position[i]<0)
			return false;
	}
	for(unsigned short i=0; i<RSA_ANGULAR_DIMENSION; i++){
		if (orientation[i]>=this->angularSize || orientation[i]<0)
			return false;
	}
	return true;

}


void PackingGenerator::testPacking(const Packing &packing, double maxTime){

	int loop = 1000*_OMP_MAXTHREADS;

	std::cout.precision(std::numeric_limits< double >::max_digits10);

	std::cout << "[" << this->seed << " PackingGenerator::testPacking] using up to " << _OMP_MAXTHREADS;
	std::cout << " concurrent treads" << std::endl;

	RND rnd(this->seed);
	RSAShape *s = ShapeFactory::createShape(&rnd);
	double dt = s->getVolume(RSA_SPATIAL_DIMENSION) / this->surface->getArea();
	delete s;

	this->packing = packing;
	for(const RSAShape *s : packing)
		this->surface->add(s);
	std::cout << "[" << this->seed << " PackingGenerator::testPacking] " << packing.size() << " shapes restored" << std::endl;



	double t = 0;

    ThreadLocalRND threadRND(rnd);

	while (t<maxTime) {

		std::cout << "\r" << "[" << this->seed << " PackingGenerator::testPacking] t=" << std::setprecision(4) << t/maxTime << " choosing " << loop << " shapes..." << std::flush;

		_OMP_PARALLEL_FOR
		for(int i = 0; i<loop; i++){
			RSAShape *sVirtual = ShapeFactory::createShape(threadRND.get());
			RSAVector pos;
			RSAOrientation angle{};
			do{
				for(unsigned short k=0; k<RSA_SPATIAL_DIMENSION; k++)
					pos[k] = threadRND.get()->nextValue()*this->spatialSize;
				for(unsigned short k=0; k<RSA_ANGULAR_DIMENSION; k++)
					angle[k] = threadRND.get()->nextValue()*this->angularSize;
			}while(!this->isInside(pos, angle));
			// setting shape position and orientation
			sVirtual->translate(pos);
			sVirtual->rotate(angle);
			// checking if shape overlaps with any shape in the packing
			if (this->surface->check(sVirtual)== nullptr){
				_OMP_CRITICAL(stdout)
				{
					std::cout << std::endl << "\t non overlapping shape found " << std::setprecision(10) << sVirtual->toString() << std::endl << std::flush;
					std::cout << "\t povray: " << std::endl << std::setprecision(10) << sVirtual->toPovray() << std::endl << std::flush;
				}
				RSAVector position = sVirtual->getPosition();
				RSAOrientation orientation = sVirtual->getOrientation();
				double delta = 0.0001;
				for(unsigned short j = 0; j< RSA_ANGULAR_DIMENSION; j++)
					orientation[j] -= 0.5*delta;

				const RSAShape *sCovers = nullptr;
				std::vector<const RSAShape*> vNeighbours;
				this->surface->getNeighbours(&vNeighbours, position);
				for(const RSAShape *sTmp : vNeighbours){
					if (sTmp->voxelInside(this->surface, position, orientation, 0.0001, delta)){
						sCovers = sTmp;
						break;
					}
				}
				if (sCovers!= nullptr)
				_OMP_CRITICAL(stdout)
				std::cout << "\t in exclusion zone of " << sCovers->toString() << std::endl;
			}
			delete sVirtual;
		} // parallel for

		t += dt * loop;
	} // while

	std::cout << "[" << seed << " PackingGenerator::testPacking] finished after time " << t << std::endl;
}


void PackingGenerator::createPacking() {

	std::cout.precision(std::numeric_limits< double >::max_digits10);
	std::cout << "[" << this->seed << " PackingGenerator::createPacking] using up to " << _OMP_MAXTHREADS;
	std::cout << " concurrent treads" << std::endl;

	std::size_t checkedAgain = 0;
	std::size_t added = 0;
	std::size_t missCounter = 0;
	unsigned short depthAnalyze = 0;

	RND rnd(this->seed);
	RSAShape *s = ShapeFactory::createShape(&rnd);
	double dt = s->getVolume(RSA_SPATIAL_DIMENSION) / this->surface->getArea();
	delete s;

	int l = 0;
	double t = 0, factor = 1;
	std::size_t tmpSplit = this->params.split, oldTmpSplit = tmpSplit;
//	int snapshotCounter = 0;

    ThreadLocalRND threadRND(rnd);

	const RSAShape **sOverlapped = new const RSAShape*[tmpSplit];
	RSAShape **sVirtual = new RSAShape*[tmpSplit];
	Voxel **aVoxels = new Voxel *[tmpSplit];

	while (!this->generationCompleted(missCounter, t)) {
		std::cout << "[" << this->seed << " PackingGenerator::createPacking] choosing " << tmpSplit << " shapes..." << std::flush;
		factor = this->getFactor();
		factor = (factor < 1.0)?1.0:factor;

		_OMP_PARALLEL_FOR
		for(std::size_t i = 0; i<tmpSplit; i++){
			sVirtual[i] = ShapeFactory::createShape(threadRND.get());
			RSAVector pos;
			RSAOrientation angle{};
			do{
				aVoxels[i] = this->voxels->getRandomVoxel(threadRND.get());
				this->voxels->getRandomPositionAndOrientation(&pos, &angle, aVoxels[i], threadRND.get());
			}while(!this->isInside(pos, angle));
			// setting shape position and orientation
			sVirtual[i]->translate(pos);
			sVirtual[i]->rotate(angle);
			// checking if shape overlaps with any shape in the packing
			sOverlapped[i] = this->surface->check(sVirtual[i]);
		} // parallel for


		std::cout << " processing shapes..." << std::flush;

		checkedAgain = 0;
		added = 0;

		// sequentially processing potentially non overlaping shapes
		for(std::size_t i=0; i<tmpSplit; i++){

			t += factor * dt;

			// if there were no intersecting particles in the packing
			if (sOverlapped[i]==nullptr){
				checkedAgain++;
				// checking overlapping again
				sOverlapped[i] = this->surface->check(sVirtual[i]);

				// if there is still no overlap sVirtula is added to the packing and corresponding voxel is removed
				if(sOverlapped[i]== nullptr){ // if non overlapping

					added++;
					l++;
					sVirtual[i]->no = l;
					sVirtual[i]->time = t;


					if (this->params.modifiedRSA){
						this->modifiedRSA(sVirtual[i], aVoxels[i]);
					}

					// consistency check
					if (aVoxels[i]!=this->voxels->getVoxel(aVoxels[i]->getPosition(), aVoxels[i]->getOrientation())){
						Voxel *v = this->voxels->getVoxel(aVoxels[i]->getPosition(), aVoxels[i]->getOrientation());
						std::cout << std::endl << "Problem: PackingGenerator - inconsistent voxels positions: [" << aVoxels[i]->toString() << "], [" << v->toString() << "]" << std::endl;
						std::cout << "size: " << this->voxels->getSpatialVoxelSize() << ", angular size: " << this->voxels->getAngularVoxelSize() << std::endl;
						std::cout << "shape: " << sVirtual[i]->toString() << std::endl;

					}

					this->surface->add(sVirtual[i]);
					this->packing.addShape(sVirtual[i]);
					this->voxels->removeTopLevelVoxel(aVoxels[i]);
				}  //second overlapping check
			} // first overlapping check
			if (sOverlapped[i]!= nullptr){ // removing overlapped virtual shapes
				delete sVirtual[i];
			}
		}  // for

		std::cout << "done, double checked: " << checkedAgain << " added: " << added << ", time: " << t << ", shapes: " << l << std::endl << std::flush;

		//whether splitting voxels
		if (added == 0) { // v.getMissCounter() % iSplit == 0){ //
			missCounter += tmpSplit;
			size_t v0 = this->voxels->getLength(), v1 = v0;

			std::cout << "[" << this->seed << " PackingGenerator::createPacking] splitting " << v0 << " voxels ";
			std::cout.flush();
//						this->toPovray("snapshot_before_" + std::to_string(snapshotCounter++) + ".pov");

			double voxelRatio = voxels->splitVoxels(this->params.minDx, this->params.maxVoxels, this->surface->getNeighbourGrid(), this->surface);
			bool b = (voxelRatio>0);
			if (b){
				v1 = this->voxels->getLength();
				v0 = (size_t)(v1/voxelRatio);
//				this->toPovray("snapshot_after_" + std::to_string(snapshotCounter++) + ".pov");
				std::cout << " done. " << this->packing.size() << " shapes, " << v1 << " voxels, new voxel size: " << voxels->getSpatialVoxelSize() << ", angular size: " << this->voxels->getAngularVoxelSize() << ", factor: " << this->getFactor() << std::endl;
				missCounter = 0;
			}else if(RSAShape::getSupportsSaturation() || rnd.nextValue() < 0.1){
				std::cout << " skipped, analyzing " << this->voxels->getLength() << " voxels, depth = " << depthAnalyze << " " << std::flush;
				this->voxels->analyzeVoxels(this->surface, this->surface->getNeighbourGrid(), depthAnalyze);
				std::cout << " done: " << this->voxels->getLength() << " voxels remained, factor = " << this->getFactor() << std::endl << std::flush;
				tmpSplit = (int)(1.1 * tmpSplit);
				v1 = this->voxels->getLength();
			}else{
				std::cout << "skipped" << std::endl << std::flush;
				depthAnalyze = 1;
			}
			// if number of voxels has changed v0!=1 is for the first split
			if (v1!=v0){
				tmpSplit *= ((double)v1 / (double)v0);
			}else{
				tmpSplit = (int)1.1*tmpSplit + _OMP_MAXTHREADS;
			}

			if (tmpSplit > std::max(this->params.maxVoxels/20, 10*this->params.split))
				tmpSplit = std::max(this->params.maxVoxels/20, 10*this->params.split);
			if(v1<v0 && voxels->getLength()<0.001*this->params.maxVoxels && tmpSplit > 10ul*_OMP_MAXTHREADS)
				tmpSplit /= 10.0;
			if(tmpSplit < 10ul*_OMP_MAXTHREADS)
				tmpSplit = 10ul*_OMP_MAXTHREADS;

			if (!b && (double)(v0-v1)/(double)v0 < 0.1){ // not much voxels removed
				depthAnalyze++;
			}else{
				if (depthAnalyze>0)
					depthAnalyze--;
			}

			if(tmpSplit != oldTmpSplit){

				delete[] sOverlapped;
				delete[] sVirtual;
				delete[] aVoxels;

				sOverlapped = new const RSAShape*[tmpSplit];
				sVirtual = new RSAShape*[tmpSplit];
				aVoxels = new Voxel *[tmpSplit];

				oldTmpSplit = tmpSplit;
			}
//			this code is for debugging purposes in case of problems with subsequent stages of packing generation

//			this->printRemainingVoxels("voxels_" + std::to_string(this->voxels->getVoxelSize()));
//			this->toWolfram("test_" + std::to_string(this->voxels->getVoxelSize()) + ".nb");
//			this->toPovray("test_" + std::to_string(this->voxels->getVoxelSize()) + ".pov");
//			this->toPovray(this->packing, this->params.surfaceSize, nullptr, "test_" + std::to_string(this->voxels->getVoxelSize()) + ".pov");
//			std::string filename = "snapshot_" + std::to_string(this->packing.size()) + "_" + std::to_string(this->voxels->getLength()) + ".dbg";
//			std::ofstream file(filename, std::ios::binary);
//			this->store(file);
//			file.close();

		}else{
			missCounter = 0;
		}
	} // while

	delete[] sOverlapped;
	delete[] sVirtual;
	delete[] aVoxels;

	std::cout << "[" << seed << " PackingGenerator::createPacking] finished after generating " << l << " shapes" << std::endl;
}

bool PackingGenerator::generationCompleted(size_t missCounter, double t) {
    bool maxDensityAchieved = false;
    if (this->params.maxDensity < 1.0) {
        double maxParticlesVolume = this->params.sufraceVolume() * this->params.maxDensity;
        maxDensityAchieved = (this->packing.getParticlesVolume(this->params.surfaceDimension) >= maxParticlesVolume);
    }

    return this->isSaturated() ||
           t >= params.maxTime ||
           missCounter >= params.maxTriesWithoutSuccess ||
           maxDensityAchieved;
}

void PackingGenerator::run(){
	this->createPacking();
}


const Packing &PackingGenerator::getPacking(){
	return this->packing;
}


void PackingGenerator::toPovray(Packing packing, double size, VoxelList *voxels, bool drawPBC,
                                const std::string &filename) {
	std::ofstream file(filename);

    if (drawPBC)
        packing.expandOnPBC(size);

    file << "#include \"colors.inc\"" << std::endl;
	file << "background { color White }" << std::endl;
	file << "camera {" << std::endl;
	file << "  orthographic" << std::endl;
	file << "  location <" << size / 2 << ", " << size / 2 << ", " << (1.3 * size) << ">" << std::endl;
	file << "  look_at <" << size / 2 << ", " << size / 2 << ",  0>" << std::endl;
	file << "  right <" << size << ", 0, 0>" << std::endl;
	file << "  up <0, " << size << ", 0>" << std::endl;
	file << "}" << std::endl;
	file << "light_source { < 1000.0, 1000.0, 1000.0> color White shadowless parallel point_at <" << size / 2 << ", " << size / 2 << ",  0>}" << std::endl;
	file << "#declare layer=union{" << std::endl;

//	file << "  polygon {5, <0.0, 0.0, 0.0>, <0.0, " << size << ", 0.0>, <" << size << ", " << size << ", 0.0>, <" << size << ", 0.0, 0.0>, <0.0, 0.0, 0.0>  texture { finish { ambient 1 diffuse 0 } pigment { color Gray} } }" << std::endl;
//	file << "  text { ttf \"timrom.ttf\" \"0\" 1, 0 pigment { color Black } scale 1.0 translate < 0, 0, 0.0002> }" << std::endl;

	for (const RSAShape *s : packing) {
//		file << "  text { ttf \"timrom.ttf\" \"" << s->no << "\" 1, 0 pigment { color White } scale 0.2 translate < " << da[0] << ", " << da[1] << ", 0.01> }" << std::endl;
//		RSAVector pos = s->getPosition();
//		if (pos[2]>0 && pos[2]<3.0)
		file << s->toPovray();
	}


	file << "}" << std::endl;

	if(voxels!=nullptr){
		file << "#declare voxels=union{" << std::endl;
		file << voxels->toPovray() << std::endl;
		file << "}" << std::endl;
	}
	file << "#declare result=union{" << std::endl;
	file << "  object { layer }" << std::endl;

	if (voxels!=nullptr)
		file << "  object { voxels }" << std::endl;
	file << "}" << std::endl;
	file << "object{ result	rotate x*360*clock }" << std::endl;

	file.close();
}


void PackingGenerator::toPovray(const std::string &filename){
    PackingGenerator::toPovray(this->packing, this->params.surfaceSize, this->voxels, false, filename);
}


void PackingGenerator::toWolfram(const RSAVector &da, const std::string &filename){

	std::vector<const RSAShape*> vShapes;
	this->surface->getNeighbours(&vShapes, da);

	std::vector<Voxel *> vVoxels;
	this->voxels->getNeighbours(&vVoxels, da);

	std::ofstream file(filename);
	file << "Graphics[{Red";

	for (const RSAShape *s : vShapes) {
		file << ", " << std::endl << s->toWolfram();
	}
	if (!vVoxels.empty()){
		file << ", Black, " << std::endl;
		for(Voxel *v: vVoxels){
			file << v->toWolfram(this->voxels->getSpatialVoxelSize(), this->voxels->getAngularVoxelSize());
		}
	}
	file << std::endl << "}]" << std::endl;
	file.close();
}


void PackingGenerator::printRemainingVoxels(const std::string &prefix){
	if (this->voxels->getLength()>20)
		return;
	for(size_t i=0; i<this->voxels->getLength(); i++){
		std::string filename(prefix + "_" + std::to_string(i) + ".nb");
		this->toWolfram(this->voxels->getVoxel(i)->getPosition(), filename);
	}
}


void PackingGenerator::toWolfram(Packing packing, double size, VoxelList *voxels, bool isPeriodicImage,
                                 const std::string &filename) {
	std::ofstream file(filename);

	if (isPeriodicImage)
	    packing.expandOnPBC(size);

    #if RSA_SPATIAL_DIMENSION == 1
		file << "Graphics[{Red";
    #elif RSA_SPATIAL_DIMENSION == 2
	    file << "Graphics[{Red";
    #elif RSA_SPATIAL_DIMENSION == 3
        file << "Graphics3D[{Red";
    #else
        die("Only 1D, 2D and 3D shapes are supported");
    #endif

	for (const RSAShape *s : packing) {
		file << ", " << std::endl << s->toWolfram();
	}

	if(voxels!=nullptr){
		file << ", Black, " << std::endl << voxels->toWolfram();
	}

	file << "}";
    if (isPeriodicImage)
        file << ", " << std::endl << "PlotRange->{{0," << size << "},{0," << size << "}}]" << std::endl;
    else
        file << "]" << std::endl;

    file.close();
}

void PackingGenerator::toWolfram(const std::string &filename){
    PackingGenerator::toWolfram(this->packing, this->params.surfaceSize, this->voxels, false, filename);
}

void PackingGenerator::store(std::ostream &f) const{
	unsigned short sd = RSA_SPATIAL_DIMENSION;
	unsigned short ad = RSA_ANGULAR_DIMENSION;
	f.write((char *)(&sd), sizeof(unsigned char));
	if (ad>0)
		f.write((char *)(&ad), sizeof(unsigned char));
	std::size_t size = this->packing.size();
	f.write((char *)(&size), sizeof(int));

	for(const RSAShape *s: this->packing){
		s->store(f);
	}
	this->voxels->store(f);
}


void PackingGenerator::restore(std::istream &f){
	unsigned char sd = RSA_SPATIAL_DIMENSION;
	unsigned char ad = RSA_ANGULAR_DIMENSION;

	f.read((char *)(&sd), sizeof(unsigned char));
	if (ad > 0)
		f.read((char *)(&ad), sizeof(unsigned char));

	if (sd!=RSA_SPATIAL_DIMENSION || ad!=RSA_ANGULAR_DIMENSION){
		std::cout << "[ERROR] cannot restore PackingGenerator: incompatible dimensions: read " << f.gcount() << " bytes." << std::endl;
		return;
	}
	int size;
	RND rnd;
	f.read((char *)(&size), sizeof(int));
	this->packing.clear();
	this->surface->clear();
	for(int i=0; i<size; i++){
		RSAShape *s = ShapeFactory::createShape(&rnd);
		s->restore(f);
		this->surface->add(s);
		this->packing.addShape(s);
	}
	this->voxels->restore(f);
}

std::vector<std::string> PackingGenerator::findPackingsInDir(const std::string &dirName) {
    std::string prefix = "packing";
    std::string suffix = ".bin";

    DIR *dir = opendir(dirName.c_str());
    dirent *de;
    std::vector<std::string> filenames;
    while ((de = readdir(dir)) != nullptr) {
        std::string filename = de->d_name;
        if (startsWith(filename, prefix) && endsWith(filename, suffix))
            filenames.push_back(dirName + "/" + filename);
    }
    (void) closedir(dir);
    return filenames;
}

std::string PackingGenerator::getPackingFilename() const {
    return this->params.getPackingSignature() + "_" + std::to_string(this->seed) + ".bin";
}
