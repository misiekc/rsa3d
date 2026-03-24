//
// Created by Michal Ciesla on 05/03/2026.
//

#include "SubVoxelList.h"

#include "utils/OMPMacros.h"



SubVoxelList::SubVoxelList(const VoxelList &vl, size_t beginIndex, size_t endIndex) : VoxelList(vl, beginIndex, endIndex) {
	this->parents = new size_t[this->length];
	this->nulls = new bool[this->length];
	this->initialLength = this->length;
	_OMP_PARALLEL_FOR
	for (size_t i=0; i<this->length; i++) {
		if (this->voxels[i]==nullptr) {
			parents[i] = this->length+1;
			this->nulls[i] = true;
		}else {
			parents[i] = i;
			this->nulls[i] = false;
		}
	}
}

SubVoxelList::SubVoxelList(const VoxelList &vl, size_t topIndex) : VoxelList(vl, 0, vl.getLength()) {
	this->parents = new size_t[this->length];
	this->nulls = new bool[this->length];
	this->initialLength = this->length;

	_OMP_PARALLEL_FOR
	for (size_t i=0; i<this->length; i++) {
		Voxel *v = this->voxels[i];
		if (this->voxels[i]==nullptr) {
			parents[i] = this->length+1;
			this->nulls[i] = true;
		}else if (this->VoxelList::getIndexOfTopLevelVoxel(v->getPosition())!=topIndex) {
			delete v;
			this->voxels[i] = nullptr;
			parents[i] = this->length+1;
			this->nulls[i] = true;
		}else {
			parents[i] = i;
			this->nulls[i] = false;
		}
	}
	this->restoreStructure();
//	std::cout << "      " << this->length << ", " << this->initialLength << " (" << this->countActiveTopLevelVoxels() << ")         ";
}

SubVoxelList::~SubVoxelList() {
	delete[] this->parents;
	delete[] this->nulls;
}

unsigned short SubVoxelList::splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc, bool verbose){
	// number of created voxels
	size_t newVoxelsCounter = 0;

	size_t voxelsFactor = 1;
	if (this->spatialVoxelSize>0)
		voxelsFactor *= (size_t)round( pow(2, this->surfaceDimension));
	RSAOrientation newAngularVoxelSize = this->angularVoxelSize;
	for (unsigned short int i=0; i<RSA_ANGULAR_DIMENSION; i++) {
		if (this->angularVoxelSize[i]>0)
			voxelsFactor *= 2;
		newAngularVoxelSize[i] /= 2;
	}

	size_t newListSize = voxelsFactor*(this->length);

	Voxel** newList = new Voxel*[ newListSize ](); // automatic nullptr initialization
	size_t *newParents = new size_t[ newListSize ];
	_OMP_PARALLEL_FOR
	for (size_t i=0; i<newListSize; i++) {
		newParents[i] = this->initialLength+1;
	}
	// temporary metrices for voxels after divisions. One separate matrix for each thread
	Voxel ***aVoxels = new Voxel**[_OMP_MAXTHREADS];
	for(int i=0; i<_OMP_MAXTHREADS; i++){
		aVoxels[i] = new Voxel*[ voxelsFactor ];
	}
	size_t dotEvery = (this->length/100)+1;
	_OMP_PARALLEL_FOR
	for(size_t i=0; i<this->length; i++){
		// voxel is tested if it should remain active and if so it is divided
		if (!this->analyzeVoxel(this->voxels[i], nl, bc, this->spatialVoxelSize, this->angularVoxelSize, 0)){
			// if too much new voxels there is no point in further splitting
			if (newVoxelsCounter <= maxVoxels){
				// preparing array of new voxels after division of this->voxels[i] (aVoxels)
				this->splitVoxel(this->voxels[i], this->spatialVoxelSize/2.0, newAngularVoxelSize, aVoxels[_OMP_THREAD_ID]);
				// analyzing new voxels, only the non covered ones will be added to the new array (newList)
				for(size_t j=0; j<voxelsFactor; j++){
					Voxel *v = aVoxels[_OMP_THREAD_ID][j];
					if( nl==nullptr || bc==nullptr || !this->analyzeVoxel(v, nl, bc, this->spatialVoxelSize/2.0, newAngularVoxelSize, 0) ){
						if(this->voxels[i]->depth > 0){
							v->depth = this->voxels[i]->depth-1;
						}
						newList[i*voxelsFactor + j] = v;
						newParents[i*voxelsFactor + j] = this->parents[i];
						_OMP_ATOMIC
						newVoxelsCounter++;
					}else{
						// covered voxels are removed
						delete aVoxels[_OMP_THREAD_ID][j];
					}
				} //for
			} // if
		}else{
			// original covered voxels are cleaned
			delete this->voxels[i];
			this->voxels[i] = nullptr;
			this->parents[i] = this->initialLength+1;
		}
		if (verbose && i%dotEvery == 0){ std::cout << "." << std::flush; }
	}
	// delete temporary thread matrices. Covered voxels have been already removed
	for(int i=0; i<_OMP_MAXTHREADS; i++){
		delete[] aVoxels[i];
	}
	delete[] aVoxels;

	if (newVoxelsCounter > maxVoxels){
		// too much voxels. Voxel splitting cancelled - using oryginal list instead of the new one
		if (verbose)
			std::cout << " too many new voxels (" << newVoxelsCounter << ">" << maxVoxels <<"): - cancel splitting," << std::flush;
		_OMP_PARALLEL_FOR
		for(size_t i=0; i<newListSize; i++){
			if (newList[i]!=nullptr)
				delete newList[i];
		}
		delete[] newList;
		delete[] newParents;

		if (verbose)
			std::cout << " compacting" << std::flush;
		this->restoreStructure();
		return VoxelList::NO_SPLIT_DUE_TO_VOXELS_LIMIT;

	}else{
		// clearing the old list and processing the new one
		_OMP_PARALLEL_FOR
		for(size_t i=0; i<this->length; i++){
			if(this->voxels[i]!=nullptr) {
				delete this->voxels[i];
				this->parents[i] = this->initialLength+1;
			}
		}
		delete[] this->voxels;
		delete[] this->parents;
		this->voxels = newList;
		this->parents = newParents;
		this->length = newListSize;

		this->spatialVoxelSize = (this->spatialVoxelSize/2.0)*this->dxFactor;

		for (unsigned short int i=0; i<RSA_ANGULAR_DIMENSION; i++) {
			this->angularVoxelSize[i] = (this->angularVoxelSize[i]/2.0)*this->dxFactor;
		}

		if (verbose)
			std::cout << " compacting" << std::flush;
		this->restoreStructure();
		return VoxelList::NORMAL_SPLIT;

	}
}

void SubVoxelList::moveVoxelInList(size_t from, size_t to){
	Assert(this->voxels[to] == nullptr);
	Assert(this->parents[to] == this->initialLength+1);
	Assert(this->voxels[from] != nullptr);
	Assert(this->parents[from] < this->initialLength);
	this->voxels[to] = this->voxels[from];
	this->parents[to] = this->parents[from];
	this->voxels[from] = nullptr;
	this->parents[from] = this->initialLength+1;
}

std::vector<size_t> SubVoxelList::getIndicesOfRemovedVoxels() const {

//	std::cout << std::endl << "     " << this->length << ", " << this->initialLength << " (" << this->countActiveTopLevelVoxels() <<")         ";

	bool *bExists = new bool[this->initialLength]();

	_OMP_PARALLEL_FOR
	for(size_t i=0; i<this->length; i++) {
		if (this->parents[i]<this->initialLength)
			bExists[this->parents[i]] = true;
	}

	std::vector<size_t> indices;
	for(size_t i=0; i<this->initialLength; i++) {
		// removed are voxels which aren't referred and wasn't null
		if(!bExists[i] && !this->nulls[i]) {
			indices.push_back(i);
		}
	}
	delete [] bExists;
	return indices;
}
