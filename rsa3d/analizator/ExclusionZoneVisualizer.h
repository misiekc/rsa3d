/*
 * ExclusionZoneVisualizer.h
 *
 *  Created on: 06.08.2018
 *      Author: ciesla
 */

#ifndef ANALIZATOR_EXCLUSIONZONEVISUALIZER_H_
#define ANALIZATOR_EXCLUSIONZONEVISUALIZER_H_

#include "../shape/Shape.h"
#include "../Parameters.h"
#include "../Packing.h"
#include "../VoxelList.h"
#include "../Surface.h"
#include "../BoundaryConditions.h"
#include <string>

class ExclusionZoneVisualizer {
private:
	const Parameters *params;
	Packing packing;
	Surface *surface;
	double voxelSize;

	bool voxelInsideShapeExclusionZone(Voxel *v, const RSAShape *shape, RSAShape *trialShape);
	bool voxelInsideExclusionZone(Voxel *v, RSAShape *trialShape);
	void normalizeShapePosition(RSAShape &shape);
	void normalizeShapeOrientation(RSAShape &shape);

public:
	ExclusionZoneVisualizer(const Parameters &params);
	void initPacking(std::string &packingFile);
	virtual ~ExclusionZoneVisualizer();
	void toPovray(const std::string &filename, RSAShape *trialShape);

	static void main(const Parameters &params, std::string &packingFile);
};

#endif /* ANALIZATOR_EXCLUSIONZONEVISUALIZER_H_ */
