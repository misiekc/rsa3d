/*
 * VariableSizeVoxelList.h
 *
 *  Created on: Nov 4, 2019
 *      Author: ciesla
 */

#ifndef VARIABLESIZEVOXELLIST_H_
#define VARIABLESIZEVOXELLIST_H_

#include "VoxelList.h"

class VariableSizeVoxelList: public VoxelList {
private:
	unsigned short *voxelsDivisionCounters;

	std::uniform_real_distribution<double> *u01Distribution;


	double *voxelMap;

	double getSpatialVoxelSize(size_t i) const;
	double getAngularVoxelSize(size_t i) const;
	double getVoxelVolume(size_t i) const;
	void createVoxelMap();

protected:

	virtual void allocateVoxels(size_t size) override;

	void moveVoxelInList(size_t from, size_t to) override;

public:
	VariableSizeVoxelList(int dim, double packingSpatialSize, double requestedSpatialVoxelSize, double shapeAngularRange, double requestedAngularVoxelSize);
	virtual ~VariableSizeVoxelList();

	void getRandomEntry(RSAVector *position, RSAOrientation *orientation, Voxel **v, RND *rnd) override;
	void getRandomPositionAndOrientation(RSAVector *position, RSAOrientation *orientation, Voxel *v, RND *rnd) override;

	unsigned short splitVoxels(double minDx, size_t maxVoxels, NeighbourGrid<const RSAShape> *nl, RSABoundaryConditions *bc) override;
	double getVoxelsVolume() const override;

	Voxel* getVoxel(const RSAVector &pos, const RSAOrientation &angle) const override;


};

#endif /* VARIABLESIZEVOXELLIST_H_ */
