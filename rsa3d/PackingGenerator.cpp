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

#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

#include "PackingGenerator.h"
#include "boundary_conditions/PeriodicBC.h"
#include "boundary_conditions/FreeBC.h"
#include "shape/ShapeFactory.h"
#include "shape/ConvexShape.h"
#include "shape/shapes/Sphere.h"
#include "ThreadLocalRND.h"
#include "utils/OMPMacros.h"
#include "surface_functions/FlatSurfaceFunction.h"
#include "surface_functions/SineSurfaceFunction.h"
#include "surface_functions/VirtualSineSurfaceFunction.h"
#include "CurvedSurface.h"
#include "CurvedSurfaceVoxelList.h"


double PackingGenerator::FACTOR_LIMIT = 5.0;


PackingGenerator::PackingGenerator(int seed, std::size_t collector, const Parameters *params)
        : seed{seed}, collector{collector}, params{*params}
{
    this->spatialSize = this->params.surfaceSize;
    this->angularSize = RSAShape::getAngularVoxelSize();
	for (unsigned short int i=0; i<RSA_ANGULAR_DIMENSION; i++) {
		if (this->params.angularVoxelRange[i] < this->angularSize[i])
			this->angularSize[i] = this->params.angularVoxelRange[i];
		if (this->params.requestedAngularVoxelSize[i] > this->angularSize[i])
			this->params.requestedAngularVoxelSize[i] = this->angularSize[i];
	}
    double gridSize = RSAShape::getNeighbourListCellSize();
    if (gridSize < this->params.thresholdDistance)
        gridSize = this->params.thresholdDistance;

    std::unique_ptr<RSABoundaryConditions> bc;
    if (this->params.boundaryConditions == "free")
        bc = std::make_unique<RSAFreeBC>();
    else
        bc = std::make_unique<RSAPeriodicBC>(this->params.surfaceSize);

    if (this->params.surfaceFunction.empty()) {
        this->surface = new Surface(this->params.surfaceDimension, this->params.surfaceSize, gridSize,
                                    RSAShape::getVoxelSpatialSize(), std::move(bc));

        this->voxels = ShapeFactory::createVoxelList(params->particleType, params->particleAttributes, this->params.surfaceDimension, this->spatialSize,
                                                     RSAShape::getVoxelSpatialSize(), this->angularSize,
                                                     this->params.requestedAngularVoxelSize);
    } else {
        Assert(this->params.surfaceDimension == RSA_SPATIAL_DIMENSION);

        std::istringstream surfaceFunctionStream(this->params.surfaceFunction);
        std::string surfaceFunctionName;
        surfaceFunctionStream >> surfaceFunctionName;
        Assert(surfaceFunctionStream);

        std::unique_ptr<SurfaceFunction> surfaceFunction;
        if (surfaceFunctionName == "FlatSurfaceFunction") {
            surfaceFunction = std::make_unique<FlatSurfaceFunction>();
        } else if (surfaceFunctionName == "SineSurfaceFunction") {
            double A{};
            std::size_t periods{};
            surfaceFunctionStream >> A >> periods;
            ValidateMsg(surfaceFunctionStream, "Malformed SineSurfaceFunctionParameters, usage: [amplitude] "
                                               "[number of periods]");
            Validate(periods > 0);
            double k = 2*M_PI*periods/this->params.surfaceSize;
            surfaceFunction = std::make_unique<SineSurfaceFunction>(A, k);
        }  else if (surfaceFunctionName == "VirtualSineSurfaceFunction") {
            double A{};
            std::size_t periods{};
            surfaceFunctionStream >> A >> periods;
            ValidateMsg(surfaceFunctionStream, "Malformed VirtualSineSurfaceFunctionParameters, usage: "
                                               "[amplitude] [number of periods]");
            Validate(periods > 0);
            double k = 2*M_PI*periods/this->params.surfaceSize;

            surfaceFunction = std::make_unique<VirtualSineSurfaceFunction>(
                A, k, Sphere<RSA_SPATIAL_DIMENSION>::getRadius()
            );
        } else {
            die("Unknown surface function: " + surfaceFunctionName);
        }

        auto curvedSurface = new CurvedSurface(this->params.surfaceDimension, this->params.surfaceSize, gridSize,
                                               RSAShape::getVoxelSpatialSize(), std::move(surfaceFunction),
                                               std::move(bc));
        this->surface = curvedSurface;
        this->voxels = new CurvedSurfaceVoxelList(this->spatialSize, RSAShape::getVoxelSpatialSize(), this->angularSize,
                                                  this->params.requestedAngularVoxelSize, curvedSurface);
    }

    if (params->maxVoxels == 0)
        this->voxels->disable();
}


PackingGenerator::~PackingGenerator() {
	delete this->voxels;
	delete this->surface;
}


bool PackingGenerator::isSaturated() {
	return (this->voxels->getLength() == 0);
}


double PackingGenerator::getFactor() {
	double factor = this->surface->getArea() / this->voxels->getVoxelsVolume();
	return (factor < 1.0)?1.0:factor;
}

void PackingGenerator::modifiedRSA(RSAShape *s, Voxel *v){

	/*const RSAShape *sn = this->surface->getClosestNeighbour(s->getPosition());
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
	}*/
}


bool PackingGenerator::isInside(const RSAVector &position, const RSAOrientation &orientation){
	for(unsigned short i=0; i<RSA_SPATIAL_DIMENSION; i++){
		if (position[i]>=this->spatialSize || position[i]<0)
			return false;
	}
	for(unsigned short i=0; i<RSA_ANGULAR_DIMENSION; i++){
		if (orientation[i]>=this->angularSize[i] || orientation[i]<0)
			return false;
	}
	return true;

}


void PackingGenerator::testPacking(const Packing &packing, double maxTime){

	int loop = 1000*_OMP_MAXTHREADS;

	std::cout.precision(std::numeric_limits< double >::max_digits10);

	std::cout << "[" << this->collector << " PackingGenerator::testPacking] using up to " << _OMP_MAXTHREADS;
	std::cout << " concurrent treads" << std::endl;

	RND rnd(this->seed);
	RSAShape *s = ShapeFactory::createShape(&rnd);
	double dt = s->getVolume(RSA_SPATIAL_DIMENSION) / this->surface->getArea();
	delete s;

	this->packing = packing;
	for(const RSAShape *s : packing)
		this->surface->add(s);
	std::cout << "[" << this->collector << " PackingGenerator::testPacking] " << packing.size() << " shapes restored" << std::endl;



	double t = 0;

    ThreadLocalRND threadRND(rnd);

	while (t<maxTime) {

		std::cout << "\r" << "[" << this->collector << " PackingGenerator::testPacking] t=" << std::setprecision(4) << t/maxTime << " choosing " << loop << " shapes..." << std::flush;

		_OMP_PARALLEL_FOR
		for(int i = 0; i<loop; i++){
			RSAShape *sVirtual = ShapeFactory::createShape(threadRND.get());
			RSAVector pos;
			RSAOrientation angle{};
			do{
				for(unsigned short k=0; k<RSA_SPATIAL_DIMENSION; k++)
					pos[k] = threadRND.get()->nextValue()*this->spatialSize;
				for(unsigned short k=0; k<RSA_ANGULAR_DIMENSION; k++)
					angle[k] = threadRND.get()->nextValue()*this->angularSize[k];
			}while(!this->isInside(pos, angle));
			// setting shape position and orientation
			sVirtual->translate(pos);
			sVirtual->rotate(angle);
			// checking if shape overlaps with any shape in the packing
			if (this->surface->check(sVirtual)== nullptr){
				_OMP_CRITICAL(stdout)
				{
					std::cout << std::endl << "\t non overlapping shape found at (" << std::setprecision(10);
					for (unsigned char i=0; i<RSA_SPATIAL_DIMENSION; i++){
						std::cout << pos[i];
						if(i<RSA_SPATIAL_DIMENSION-1)
							std::cout << ", ";
					}
					std::cout << "), (";
					for (unsigned char i=0; i<RSA_ANGULAR_DIMENSION; i++){
						std::cout << angle[i];
						if(i<RSA_ANGULAR_DIMENSION-1)
							std::cout << ", ";
					}
					std::cout <<")" << std::endl;

					std::cout << sVirtual->toString() << std::endl << std::flush;
					std::cout << "\t povray: " << std::endl << std::setprecision(10) << sVirtual->toPovray() << std::endl << std::flush;
				}
				RSAVector position = sVirtual->getPosition();
				RSAOrientation orientation = sVirtual->getOrientation();
				RSAOrientation angularRange;
				const double delta = 0.0001;
				for(unsigned short j = 0; j< RSA_ANGULAR_DIMENSION; j++) {
					orientation[j] -= 0.5*delta;
					angularRange[j] = delta;
				}

				const RSAShape *sCovers = nullptr;
				std::vector<const RSAShape*> vNeighbours;
				this->surface->getNeighbours(&vNeighbours, position);
				for(const RSAShape *sTmp : vNeighbours){
					if (sTmp->voxelInside(this->surface->getBC(), position, orientation, delta, angularRange)){
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

	std::cout << "[" << this->collector << " PackingGenerator::testPacking] finished after time " << t << std::endl;
}


void PackingGenerator::createPacking(Packing *packing) {

	std::cout.precision(std::numeric_limits< double >::max_digits10);
	std::cout << "[" << this->collector << " PackingGenerator::createPacking] using up to " << _OMP_MAXTHREADS;
	std::cout << " concurrent treads" << std::endl;

	std::size_t checkedAgain = 0;
	std::size_t added = 0;
	std::size_t addedSinceLastSplit = 0;
	std::size_t missCounter = 0;
	unsigned short depthAnalyze = 0;

	RND rnd(this->seed);
	RSAShape *s = ShapeFactory::createShape(&rnd);
	double dt = s->getVolume(RSA_SPATIAL_DIMENSION) / this->surface->getArea();
	delete s;

	int l = 0;
	double t = 0;
	double factor = this->getFactor();
	if (packing!=nullptr){
		this->packing = *packing;
		for(const RSAShape *s : *packing)
			this->surface->add(s);
		l = this->packing.size();
	}

	std::size_t tmpSplit = this->params.split, oldTmpSplit = tmpSplit;
//	int snapshotCounter = 0;

    ThreadLocalRND threadRND(rnd);

	const RSAShape **sOverlapped = new const RSAShape*[tmpSplit];
	RSAShape **sVirtual = new RSAShape*[tmpSplit];
	Voxel **aVoxels = new Voxel *[tmpSplit];

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	while (!this->generationCompleted(missCounter, t)) {
		if (this->params.timestamp){
		    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
		    auto miliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - begin).count();
		    std::cout << "[" << miliseconds << "]";

		}
		std::cout << "[" << this->collector << " PackingGenerator::createPacking] choosing " << tmpSplit << " shapes..." << std::flush;

		_OMP_PARALLEL_FOR
		for(std::size_t i = 0; i<tmpSplit; i++){
			sVirtual[i] = ShapeFactory::createShape(threadRND.get());
			RSAVector pos;
			RSAOrientation angle{};
			do{
			    this->voxels->getRandomEntry(&pos, &angle, &(aVoxels[i]), threadRND.get());
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
						std::cout << "size: " << this->voxels->getSpatialVoxelSize() << ", angular size: ";
						RSAOrientation avSize = this->voxels->getAngularVoxelSize();
						for (unsigned short int j=0; j<RSA_ANGULAR_DIMENSION; j++) {
							std::cout << avSize[j] << ", ";
						}
						std::cout << '\b' << '\b' << std::endl;
						std::cout << ", shape: " << sVirtual[i]->toString() << std::endl;
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
		addedSinceLastSplit += added;
		//whether splitting voxels
		if (added == 0) { // v.getMissCounter() % iSplit == 0){ //
			missCounter += tmpSplit;
			size_t v0 = this->voxels->getLength();

			if (this->params.timestamp){
			    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
			    auto miliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - begin).count();
			    std::cout << "[" << miliseconds << "]";

			}
			std::cout << "[" << this->collector << " PackingGenerator::createPacking] shapes added since last split " << addedSinceLastSplit << std::endl;
			std::cout << "[" << this->collector << " PackingGenerator::createPacking] splitting " << v0 << " voxels ";
			std::cout.flush();
//						this->toPovray("snapshot_before_" + std::to_string(snapshotCounter++) + ".pov");

			unsigned short status = voxels->splitVoxels(this->params.minDx, this->params.maxVoxels, this->surface->getNeighbourGrid(), this->surface->getBC());
			double oldFactor = factor;
			factor = this->getFactor();

			if (status == VoxelList::NORMAL_SPLIT || status == VoxelList::NO_SPLIT_BUT_INITIALIZED){
//				this->toPovray("snapshot_after_" + std::to_string(snapshotCounter++) + ".pov");
				std::cout.precision(5);
				std::cout << " done. " << this->voxels->getLength() << " (" << this->voxels->countActiveTopLevelVoxels() << ") voxels, new voxel size: " << voxels->getSpatialVoxelSize() << ", angular size: ";
				RSAOrientation avSize = this->voxels->getAngularVoxelSize();
				for (unsigned short int i=0; i<RSA_ANGULAR_DIMENSION; i++) {
					std:: cout << avSize[i] << " ";
				}
				std::cout << ", factor: " << factor << ", change: " << (factor/oldFactor) << std::endl;
				std::cout.precision(std::numeric_limits< double >::max_digits10);
				missCounter = 0;
				addedSinceLastSplit = 0;
			}else if (status == VoxelList::NO_SPLIT_DUE_TO_VOXELS_LIMIT){
				if(RSAShape::getSupportsSaturation() || rnd.nextValue() < 0.1){
					if (!this->params.goDeep && addedSinceLastSplit==0){
						std::cout << " skipped after generating " << l << " shapes" << std::endl;
						delete[] sOverlapped;
						delete[] sVirtual;
						delete[] aVoxels;
						return;
					}
					std::cout << " skipped, analyzing " << this->voxels->getLength() << " voxels, depth = " << depthAnalyze << " " << std::flush;
					this->voxels->analyzeVoxels(this->surface->getBC(), this->surface->getNeighbourGrid(), depthAnalyze);
					factor = this->getFactor();
					std::cout.precision(5);
					std::cout << " done: " << this->voxels->getLength() << " (" << this->voxels->countActiveTopLevelVoxels() << ") voxels remained, factor = " << factor << ", change: " << (factor/oldFactor) << std::endl << std::flush;
					std::cout.precision(std::numeric_limits< double >::max_digits10);
					tmpSplit = (int)(1.1 * tmpSplit);
				}
				addedSinceLastSplit = 0;
			}else{
				std::cout << "skipped" << std::endl << std::flush;
				depthAnalyze = 1;
				addedSinceLastSplit = 0;
			}
			size_t v1 = this->voxels->getLength();

			// determining new value of tmpSplit


			if (status == VoxelList::NO_SPLIT_BUT_INITIALIZED){
				tmpSplit = (int)(tmpSplit/factor);
			}else{
				// standard grow of tmpSplit
				tmpSplit = (int)(tmpSplit * 1.1* v1 / v0);
			}
			// increase split when number of voxels closes to its limit
			if (this->voxels->getLength() > 0.05*this->params.maxVoxels)
				tmpSplit *= 1.5;
			// increase split when number of voxels limit is reached
			if (status == VoxelList::NO_SPLIT || status == VoxelList::NO_SPLIT_DUE_TO_VOXELS_LIMIT){
				if ((double)(v0-v1)/(double)v0 < 0.1){ // not much voxels removed
					depthAnalyze++;
					tmpSplit = 2*tmpSplit;
				}
			}else{
				if (depthAnalyze>0)
					depthAnalyze--;
			}

			// additional tweaking

			if (tmpSplit > std::max(this->params.maxVoxels/5, 10*this->params.split))
				tmpSplit = std::max(this->params.maxVoxels/5, 10*this->params.split);
			if(v1<v0 && voxels->getLength()<0.001*this->params.maxVoxels && tmpSplit > 10ul*_OMP_MAXTHREADS)
				tmpSplit /= 10.0;

			tmpSplit = ((tmpSplit / _OMP_MAXTHREADS)+1) * _OMP_MAXTHREADS;

			if(tmpSplit < 10ul*_OMP_MAXTHREADS)
				tmpSplit = 10ul*_OMP_MAXTHREADS;

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

//			this->printRemainingVoxels("voxels_" + std::to_string(this->voxels->getSpatialVoxelSize()));
//			this->toWolfram("test_" + std::to_string(this->voxels->getSpatialVoxelSize()) + ".nb");
//			this->toPovray("test_" + std::to_string(this->voxels->getSpatialVoxelSize()) + ".pov");
//			this->toPovray(this->packing, this->params.surfaceSize, nullptr, false, "test_" + std::to_string(this->voxels->getSpatialVoxelSize()) + ".pov");
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

	if (this->params.timestamp){
	    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
	    auto miliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - begin).count();
	    std::cout << "[" << miliseconds << "]";

	}
	std::cout << "[" << this->collector << " PackingGenerator::createPacking] finished after generating " << l << " shapes" << std::endl;
}

bool PackingGenerator::generationCompleted(size_t missCounter, double t) {
    bool maxDensityAchieved = false;
    if (this->params.maxDensity < 1.0) {
        double maxParticlesVolume = this->params.sufraceVolume() * this->params.maxDensity;
        maxDensityAchieved = (this->packing.getParticlesVolume(this->params.packingFractionSurfaceDimension())
                              >= maxParticlesVolume);
    }

    return this->isSaturated() ||
           t >= params.maxTime ||
           missCounter >= params.maxTriesWithoutSuccess ||
           maxDensityAchieved;
}


void PackingGenerator::run(Packing *packing){

#ifdef _OPENMP
	struct timeval timeout;
	timeout.tv_sec = 0;
	timeout.tv_usec = 0;
	std::string s = ""; //default to empty string
	fd_set fds;
	FD_ZERO(&fds);
	FD_SET(STDIN_FILENO, &fds);
	int result = select(STDIN_FILENO + 1, &fds, NULL, NULL, &timeout);
	if (result != -1){
		if (FD_ISSET(STDIN_FILENO, &fds)){
			std::string s;
			std::cin >> s;
			std::string sOmpThreads = "ompThreads:";
			size_t pos = s.find(sOmpThreads);
			if (pos!=std::string::npos){
				s = s.substr(sOmpThreads.length(), s.length()-sOmpThreads.length());
				omp_set_num_threads(std::stoi(s));
			}
		}
	}
#endif

	this->createPacking(packing);
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
	file << "  location <" << size / 2 << ", " << size / 2 << ", " << (1.1*size) << ">" << std::endl;
	file << "  look_at <" << size / 2 << ", " << size / 2 << ",  0>" << std::endl;
	file << "  right <" << (1.1*size) << ", 0, 0>" << std::endl;
	file << "  up <0, " << (1.1*size) << ", 0>" << std::endl;
	file << "}" << std::endl;
	file << "light_source { < 1000.0, 1000.0, 1000.0> color White shadowless parallel point_at <" << size / 2 << ", " << size / 2 << ",  0>}" << std::endl;
	file << "#declare layer=union{" << std::endl;

//	file << "  polygon {5, <0.0, 0.0, 0.0>, <0.0, " << size << ", 0.0>, <" << size << ", " << size << ", 0.0>, <" << size << ", 0.0, 0.0>, <0.0, 0.0, 0.0>  texture { finish { ambient 1 diffuse 0 } pigment { color Gray} } }" << std::endl;
//	file << "  text { ttf \"timrom.ttf\" \"0\" 1, 0 pigment { color Black } scale 1.0 translate < 0, 0, 0.0002> }" << std::endl;

	for (const RSAShape *s : packing) {
//		RSAVector pos = s->getPosition();
//		file << "  text { ttf \"timrom.ttf\" \"" << s->no << "\" 1, 0 pigment { color White } scale 0.2 translate < " << pos[0] << ", " << pos[1] << ", 0.01> }" << std::endl;
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
    return this->params.getPackingSignature() + "_" + std::to_string(this->collector) + ".bin";
}
