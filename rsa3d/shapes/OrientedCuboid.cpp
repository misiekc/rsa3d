/*
 * OrientedCuboid.cpp
 *
 *  Created on: 12.07.2017
 *      Author: ciesla
 */

#include "OrientedCuboid.h"
#include <algorithm>

bool OrientedCuboid::do2Drotation;

OrientedCuboid::OrientedCuboid(const Matrix<3, 3> & rotation) : Cuboid(rotation){
}

OrientedCuboid::~OrientedCuboid() {
}

void OrientedCuboid::initClass(const std::string &args){
    Cuboid::initClass(args);
    if (args.compare(args.length() - 4, 4, "true")==0)
    	OrientedCuboid::do2Drotation = true;
    else{
    	OrientedCuboid::do2Drotation = false;
    	neighbourListCellSize = *std::max_element(size, size + staticDimension);
    	voxelSize = 0.5*(*std::min_element(size, size + staticDimension));
    }
}

// Method creating (dynamically alocated) cuboid with random orientation.
// Used by ShapeFactory for shape generating
//----------------------------------------------------------------------------
Shape * OrientedCuboid::create(RND *rnd){
	OrientedCuboid *cuboid;
	if (OrientedCuboid::do2Drotation)
	    cuboid = new OrientedCuboid(Matrix<3, 3>::rotation(
	        0,
	        std::asin(rnd->nextValue() * 2 - 1),
	        rnd->nextValue() * 2 * M_PI));
	else
    	cuboid = new OrientedCuboid(Matrix<3, 3>::identity());

#ifdef CUBOID_DEBUG
    std::cout << "Creating OrientedCuboid:" << std::endl;
    std::cout << cuboid->orientation;
#endif

    return cuboid;
}

int OrientedCuboid::overlap(BoundaryConditions *bc, Shape *s){
    return this->pointInside(bc, s->getPosition());
}

int OrientedCuboid::pointInside(BoundaryConditions *bc, double* da){
    if (OrientedCuboid::do2Drotation){
    	return Cuboid::pointInside(bc, da);
    }else{
    	double *ta = new double[this->dimension];
   		bc->getTranslation(ta, this->position, da);
    	for(unsigned char i=0; i<this->dimension; i++){
    		if (std::fabs(this->position[i] - (da[i] + ta[i])) > size[i]){
    			delete[] ta;
    			return false;
    		}
    	}
    	delete[] ta;
    	return true;
    }
}
