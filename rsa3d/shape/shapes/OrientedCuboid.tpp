/*
 * OrientedCuboid.cpp
 *
 *  Created on: 12.07.2017
 *      Author: ciesla
 */

#include <algorithm>
#include <string>

// bool OrientedCuboid::do2Drotation;
template <unsigned short DIMENSION>
double OrientedCuboid<DIMENSION>::size[DIMENSION];

template <unsigned short DIMENSION>
void OrientedCuboid<DIMENSION>::initClass(const std::string &args){
	
    std::stringstream args_stream (args);
        
    for (unsigned short i = 0; i < DIMENSION; i++) {
        args_stream >> size[i];
        if (!args_stream)           throw std::runtime_error("OrientedCuboid::initClass: invalid or missing dimensions");
        else if (size[i] <= 0.0)    throw std::runtime_error("OrientedCuboid::initClass: non-positive size: " + std::to_string(size[i]));
    }

    // renormailze sizes to obtain unit volume
    double volume = std::accumulate(size, size + DIMENSION, 1.0, std::multiplies<double>());
    double factor = 1.0/pow(volume, 1.0/DIMENSION);
    for (unsigned short i = 0; i < DIMENSION; i++) {
        size[i] *= factor;
    }

    double diagonal = std::accumulate(size, size + DIMENSION, 0.0, [](auto squareSum, auto size) {
        return squareSum + size*size;
    });
    diagonal = std::sqrt(diagonal);
    double minSize = *std::min_element(OrientedCuboid<DIMENSION>::size, OrientedCuboid<DIMENSION>::size + DIMENSION);
    double maxSize = *std::max_element(OrientedCuboid<DIMENSION>::size, OrientedCuboid<DIMENSION>::size + DIMENSION);

    ShapeStaticInfo<DIMENSION, 0> shapeInfo;
    shapeInfo.setCircumsphereRadius(diagonal / 2);
    shapeInfo.setInsphereRadius(minSize / 2);
    shapeInfo.setNeighbourListCellSize(maxSize);
    shapeInfo.setSpatialVoxelSize(minSize);
    shapeInfo.setSupportsSaturation(true);
	shapeInfo.template setDefaultCreateShapeImpl <OrientedCuboid<DIMENSION>> ();

	Shape<DIMENSION, 0>::setShapeStaticInfo(shapeInfo);
	Shape<DIMENSION, 0>::setEarlyRejectionEnabled(false);     // it doesn't help here
}

template <unsigned short DIMENSION>
bool OrientedCuboid<DIMENSION>::overlap(BoundaryConditions<DIMENSION> *bc, const Shape<DIMENSION, 0> *s) const{
	Orientation<0> oRange;
    return this->pointInside(bc, s->getPosition(), s->getOrientation(), oRange);
}

template <unsigned short DIMENSION>
bool OrientedCuboid<DIMENSION>::pointInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &position,
										   const Orientation<0> &orientation, const Orientation<0> &orientationRange) const{
	Vector<DIMENSION> ta = bc->getTranslation(this->getPosition(), position);
	for(unsigned short i=0; i<DIMENSION; i++){
		if (std::fabs(this->getPosition()[i] - (position[i] + ta[i])) > size[i]){
			return false;
		}
	}
	return true;
}

/*
 * -1 - outside,
 * otherwise returned value contains a set of flags showin in which dimension the voxel is fully covered.
 * For example for three dimensional cuboid the value 3 means that the voxel is fully covered in x and y coordinate end partially in z
 * 0 - 2^d-1 partially inside, 2^d-1 fully inside
 */

template <unsigned short DIMENSION>
int OrientedCuboid<DIMENSION>::voxelFullyOrPartiallyInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &voxelPosition, double spatialSize) const {

	Orientation<0> orientation;
	Orientation<0> angularSize;

    switch(this->voxelInsideEarlyRejection(bc, voxelPosition, orientation, spatialSize, angularSize)) {
        case EarlyRejectionResult::TRUE:      return (1 << DIMENSION)-1;
        case EarlyRejectionResult::FALSE:     return -1;
        case EarlyRejectionResult::UNKNOWN:   break;
    }

    RSAVector position = this->getPosition();
	RSAVector ta = bc->getTranslation(position, voxelPosition);
	int iResult = 0;

    for(unsigned short i=0; i<DIMENSION; i++){
    	bool fullyInside = (voxelPosition[i] + ta[i] >= position[i] - size[i]) && (voxelPosition[i] + ta[i] + spatialSize < position[i] + size[i]);
    	if (fullyInside){
    		iResult |= (1 << i);
    	}else{
    		bool fullyOutside = (voxelPosition[i] + ta[i] > position[i] + size[i]) || (voxelPosition[i] + ta[i] + spatialSize <= position[i] - size[i]);
    		if (fullyOutside)
    			return -1;
    	}
    }
    return iResult;
}

template <unsigned short DIMENSION>
bool OrientedCuboid<DIMENSION>::voxelInside(BoundaryConditions<DIMENSION> *bc,
											const Vector<DIMENSION> &voxelPosition,
											const Orientation<0> &orientation,
                                            double spatialSize, const Orientation<0> &angularSize) const {

	return (this->voxelFullyOrPartiallyInside(bc, voxelPosition, spatialSize) == (1 << DIMENSION)-1);
}

/*
 * @brief Method used by OrientedCuboidVoxelList to chceck if voxel is overlapped by a set of neighbouring shapes
 *
 */

template <unsigned short DIMENSION>
bool OrientedCuboid<DIMENSION>::voxelInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &voxelPosition, double spatialSize, const std::vector<const RSAShape *> *shapes){
	// check over single shapes
	std::vector< std::pair<const RSAShape *, unsigned short> > closestShapes;
//	RSAVector pos;
//	RSAOrientation orientation;
	// check over shapes and select positions of shapes which exclusion zones partially overlap the voxel
	for(const RSAShape *s: *shapes){
//		if (s->voxelInside(bc, voxelPosition, orientation, spatialSize, 0)){
		int res = ((OrientedCuboid<DIMENSION> *)s)->voxelFullyOrPartiallyInside(bc, voxelPosition, spatialSize);
		if (res==(1 << DIMENSION)-1){
			return true;
		}else if (res!=-1){
			unsigned int flags = (unsigned int)res;
			unsigned int all1 = (1 << DIMENSION)-1;
			unsigned int invertedFlags = ~flags & all1;
			double ddim = std::log2( invertedFlags );
			if (ddim != round(ddim))
				continue;
			std::pair <const RSAShape *, unsigned short> p;
			p.first = s;
			p.second = (unsigned short)ddim;
			closestShapes.push_back(p);
		}
	}

	// check if pair of shapes fully overlap the voxel
	// we look for pair of such shapes which partially overlaps the voxel along the same axis, voxel are between them and their distance is lower than 2*size
	for (size_t i=0; i<closestShapes.size(); i++){
		RSAVector pi = closestShapes[i].first->getPosition() + bc->getTranslation(voxelPosition, closestShapes[i].first->getPosition());
		int k = closestShapes[i].second;
		double ipos = pi[k];
		double vpos = voxelPosition[k];
		for (size_t j=i+1; j<closestShapes.size(); j++){
			if (closestShapes[j].second!=k)
				continue; // different coordinates - not interesting
			RSAVector pj = closestShapes[j].first->getPosition() + bc->getTranslation(voxelPosition, closestShapes[j].first->getPosition());
			double jpos = pj[k];
			if ( (ipos - size[k] < vpos && vpos + spatialSize < jpos + size[k] ) && (jpos-ipos < 2*size[k]) )
				return true;
			if ( (jpos - size[k] < vpos && vpos + spatialSize < ipos + size[k] ) && (ipos-jpos < 2*size[k]) )
				return true;
		}
	}
	return false;
}


template <unsigned short DIMENSION>
double OrientedCuboid<DIMENSION>::getVolume(unsigned short dim) const{
    if (dim != DIMENSION)
        throw std::runtime_error("OrientedCuboid supports only " + std::to_string(DIMENSION) + "D packings");
    return std::accumulate(size, size + DIMENSION, 1.0, std::multiplies<double>());
}

template <unsigned short DIMENSION>
std::string OrientedCuboid<DIMENSION>::toPovray() const
{	
	std::string s = "";
	Vector<DIMENSION> position = this->getPosition();
	if (DIMENSION==2){
		double factor = 0.5;
		s += "  polygon {4, ";

		s += "< " + std::to_string(position[0] - factor*size[0]) + ", ";
		s += 		std::to_string(position[1] - factor*size[1]) + ", 0.1>, ";

		s += "< " + std::to_string(position[0] - factor*size[0]) + ", ";
		s += 		std::to_string(position[1] + factor*size[1]) + ", 0.1>, ";

		s += "< " + std::to_string(position[0] + factor*size[0]) + ", ";
		s += 		std::to_string(position[1] + factor*size[1]) + ", 0.1>, ";

		s += "< " + std::to_string(position[0] + factor*size[0]) + ", ";
		s += 		std::to_string(position[1] - factor*size[1]) + ", 0.1> ";

		s += "\n    texture { pigment { color Red } }\n  }\n";
	}
	else if (DIMENSION==3){
		s += "  box { < ";
		for(unsigned short i=0; i<DIMENSION; i++){
			s += std::to_string(-OrientedCuboid<DIMENSION>::size[i]/2);
			if (i<DIMENSION-1)
				s+= ", ";
		}
		s += ">, <";
		for(unsigned short i=0; i<DIMENSION; i++){
			s += std::to_string(OrientedCuboid<DIMENSION>::size[i]/2);
			if (i<DIMENSION-1)
				s+= ", ";
		}
		s += ">\n";
		s += "    matrix < \n    ";
		for (unsigned short i=0; i<DIMENSION; i++){
			for (unsigned short j=0; j<DIMENSION; j++){
				s += i==j?"1":"0";
				if (j<DIMENSION-1)
					s+= ", ";
			}
			s+= ",\n    ";
		}

		for(unsigned short i=0; i<DIMENSION; i++){
			s += std::to_string(position[i]);
			if (i<DIMENSION-1)
				s+= ", ";
		}
		s += "\n    >\n";
		s += "    texture { pigment { color Red } }\n  }\n";
	}
	return s;
}

template<unsigned short DIMENSION>
Shape<DIMENSION, 0> *OrientedCuboid<DIMENSION>::clone() const {
    return new OrientedCuboid(*this);
}
