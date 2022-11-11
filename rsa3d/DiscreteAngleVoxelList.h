/*
 * DiscreteAngleVoxelList.h
 *
 *  Created on: 11.11.2022
 *      Author: Michal Ciesla
 */

#ifndef DISCRETEANGLEVOXELLIST_H_
#define DISCRETEANGLEVOXELLIST_H_

#include "VoxelList.h"

class DiscreteAngleVoxelList: public VoxelList {

private:
protected:
	virtual unsigned int initVoxels(RSABoundaryConditions *bc, NeighbourGrid<const RSAShape> *nl) override;
	virtual void splitVoxel(Voxel *v, double spatialSize, double angularSize, Voxel **vRes) override;
	virtual void getRandomPositionAndOrientation(RSAVector *position, RSAOrientation *orientation, Voxel *v, RND *rnd);

public:
	std::vector<RSAOrientation> allowedOrientations;
	/**
	 * @brief Constructor for list allowing specific orientations only
	 * @param packingSpatialSize packing size to be covered by voxels
	 * @param requestedSpatialVoxelSize suggested initial size of a voxel. Initial size of a allocated voxels will not be larger than the requested one.
	 * @param orientations vector of allowed orientations.
	 */
	DiscreteAngleVoxelList(int dim, double packingSpatialSize, double requestedSpatialVoxelSize, const std::vector<Orientation<1>> &orientations);

	Voxel *getVoxel(const RSAVector &pos, const RSAOrientation &angle) const override;
	virtual double getAngularVoxelSize() const override;
	virtual double getVoxelsVolume() const override;


};


#endif /* DISCRETEANGLEVOXELLIST_H_ */
