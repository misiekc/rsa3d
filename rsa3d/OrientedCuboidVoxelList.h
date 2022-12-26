/*
 * OrientedCuboidVoxelList.h
 *
 *  Created on: 30.09.2019
 *      Author: ciesla
 */

#ifndef ORIENTEDCUBOIDVOXELLIST_H_
#define ORIENTEDCUBOIDVOXELLIST_H_

#include "VariableSizeVoxelList.h"

class OrientedCuboidVoxelList : public VariableSizeVoxelList {

	protected:
		bool isVoxelInsideExclusionZone(Voxel *v, double spatialSize, double angularSize, std::vector<const RSAShape*> *shapes, RSABoundaryConditions *bc, unsigned short depth = 0) const override ;

	public:
		OrientedCuboidVoxelList(int dim, double packingSpatialSize, double requestedSpatialVoxelSize, double shapeAngularRange, double requestedAngularVoxelSize);
};

#endif /* ORIENTEDCUBOIDVOXELLIST_H_ */
