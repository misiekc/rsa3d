/*
 * ExclusionZoneVisualizer.cpp
 *
 *  Created on: 06.08.2018
 *      Author: ciesla
 */

#include "ExclusionZoneVisualizer.h"
#include "../surfaces/NBoxPBC.h"
#include "../RND.h"
#include "../shape/ShapeFactory.h"
#include "../Voxel.h"
#include "../Vector.h"
#include <fstream>

ExclusionZoneVisualizer::ExclusionZoneVisualizer(const Parameters &params) {
	this->params = &params;
	this->surface = new NBoxPBC(this->params->surfaceDimension, this->params->surfaceSize, RSAShape::getNeighbourListCellSize(), RSAShape::getVoxelSpatialSize());
	this->voxelSize = this->params->surfaceSize/800;
}

void ExclusionZoneVisualizer::initPacking(std::string &packingFile) {
	this->packing.restore(packingFile);
	for(const RSAShape *s: packing){
		this->surface->add(s);
	}
}

ExclusionZoneVisualizer::~ExclusionZoneVisualizer() {
	delete this->surface;
}

bool ExclusionZoneVisualizer::voxelInsideShapeExclusionZone(Voxel *v, const RSAShape *shape, RSAShape *trialShape){
    RSAVector position;
    RSAVector voxelPosition = v->getPosition();
    int counterSize = 1 << RSA_SPATIAL_DIMENSION;
    bool isInside = true;
    for(int i=0; i<counterSize; i++){
        for(unsigned short j=0; j<RSA_SPATIAL_DIMENSION; j++){
            position[j] = voxelPosition[j] + RSAPositioned::offset[i][j]*this->voxelSize;
        }
        trialShape->translate(position - trialShape->getPosition());
        if(!shape->overlap(this->surface, trialShape)){
            isInside = false;
            break;
        }
    }
    return isInside;
}

bool ExclusionZoneVisualizer::voxelInsideExclusionZone(Voxel *v, RSAShape *trialShape){
	std::vector<const RSAShape *> neighbours;
	this->surface->getNeighbours(&neighbours, v->getPosition());
	for (const RSAShape *s : neighbours){
		if (this->voxelInsideShapeExclusionZone(v, s, trialShape))
			return true;
	}
	return false;
}

void ExclusionZoneVisualizer::toPovray(const std::string &filename, RSAShape *trialShape){
	std::ofstream file(filename);

	file << "#include \"colors.inc\"" << std::endl;
	file << "background { color Black }" << std::endl;
	file << "camera { orthographic location <" << this->params->surfaceSize / 2 << ", " << this->params->surfaceSize / 2 << ", " << (1.4 * this->params->surfaceSize) << "> look_at  <" << this->params->surfaceSize / 2 << ", " << this->params->surfaceSize / 2 << ",  0> }" << std::endl;
	file << "light_source { <" << this->params->surfaceSize / 2 << ", " << this->params->surfaceSize / 2 << ", " << (1.5 * this->params->surfaceSize) << "> color White shadowless parallel point_at <" << this->params->surfaceSize / 2 << ", " << this->params->surfaceSize / 2 << ",  0>}" << std::endl;
	file << "#declare layer=union{" << std::endl;

	file << "  polygon {5, <0.0, 0.0, 0.0>, <0.0, " << this->params->surfaceSize << ", 0.0>, <" << this->params->surfaceSize << ", " << this->params->surfaceSize << ", 0.0>, <" << this->params->surfaceSize << ", 0.0, 0.0>, <0.0, 0.0, 0.0>  texture { finish { ambient 1 diffuse 0 } pigment { color White} } }" << std::endl;

	RSAVector position;
/*
	position[0] = -1.0;
	position[1] = 4.5;
	trialShape->translate(position - trialShape->getPosition());
	file << trialShape->toPovray();
*/
	for (const RSAShape *s : packing) {
		file << s->toPovray();
	}


	file << "}" << std::endl;

	// voxels
	file << "#declare voxels=union{" << std::endl;

	size_t ns = (int)(this->params->surfaceSize / this->voxelSize) + 1;
	std::array<int, RSA_SPATIAL_DIMENSION> ins{};
    ins.fill(0);
	RSAOrientation orientation{};
    orientation.fill(0);

	do{
		for(unsigned char i=0; i<RSA_SPATIAL_DIMENSION; i++){
			position[i] = this->voxelSize*ins[i]; // position point to the "left bottom" corner of a voxel
		}
		Voxel v(position, orientation);
		if(this->voxelInsideExclusionZone(&v, trialShape))
			file << v.toPovray(this->voxelSize) << std::endl;
	}while(increment(ins.data(), RSA_SPATIAL_DIMENSION, ns-1));

	file << "}" << std::endl;

	file << "#declare result=union{" << std::endl;
	file << "  object { layer }" << std::endl;
	file << "  object { voxels }" << std::endl;
	file << "}" << std::endl;
	file << "object{ result	rotate x*360*clock }" << std::endl;

	file.close();
}

void ExclusionZoneVisualizer::normalizeShapePosition(RSAShape &shape){
	Vector<RSA_SPATIAL_DIMENSION, double> position{};
	position = shape.getPosition();
	for(size_t i=0; i< position.getDimension(); i++){
		position[i] = -position[i] + 0.5*this->params->surfaceSize;
	}
	shape.translate(position);
}

void ExclusionZoneVisualizer::normalizeShapeOrientation(RSAShape &shape){
	std::array<double, RSA_ANGULAR_DIMENSION> orientation{};
	orientation = shape.getOrientation();
	for(size_t i=0; i< orientation.size(); i++){
		orientation[i] = -orientation[i];
	}
	shape.rotate(orientation);
}

void ExclusionZoneVisualizer::main(const Parameters &params, std::string &packingFile, std::string &outputFile){
	ExclusionZoneVisualizer xzv(params);
	xzv.initPacking(packingFile);

	RND rnd;
  	std::array<double, RSA_ANGULAR_DIMENSION> orientation{};
/*
	RSAShape *tmpShape = ShapeFactory::createShape(&rnd);
	xzv.normalizeShapePosition(*tmpShape);
	xzv.normalizeShapeOrientation(*tmpShape);
	xzv.packing.addShape(tmpShape);
	xzv.surface->add(tmpShape);
*/
	RSAShape *trialShape = ShapeFactory::createShape(&rnd);
	xzv.normalizeShapeOrientation(*trialShape);

	orientation[0] = M_PI / 180;
	char buf[20];

	for(int i=0; i<180; i++){
		trialShape->rotate(orientation);
		std::sprintf(buf, "output_%03d.pov", i);
		xzv.toPovray(std::string(buf), trialShape);
	}
	delete trialShape;
}


