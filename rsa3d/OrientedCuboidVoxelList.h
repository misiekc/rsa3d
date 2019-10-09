/*
 * OrientedCuboidVoxelList.h
 *
 *  Created on: 30.09.2019
 *      Author: ciesla
 */

#ifndef ORIENTEDCUBOIDVOXELLIST_H_
#define ORIENTEDCUBOIDVOXELLIST_H_

#include "VoxelList.h"

class OrientedCuboidVoxelList : public VoxelList {

	protected:
		bool isVoxelInsideExclusionZone(Voxel *v, double spatialSize, double angularSize, std::vector<const RSAShape*> *shapes, RSABoundaryConditions *bc, unsigned short depth = 0) override ;

	public:
		OrientedCuboidVoxelList(int dim, double packingSpatialSize, double requestedSpatialVoxelSize, double shapeAngularRange, double requestedAngularVoxelSize);
};

#endif /* ORIENTEDCUBOIDVOXELLIST_H_ */
