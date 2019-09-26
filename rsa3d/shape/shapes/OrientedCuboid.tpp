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
    shapeInfo.setVoxelSpatialSize(minSize / 2);
    shapeInfo.setSupportsSaturation(true);
	shapeInfo.template setDefaultCreateShapeImpl <OrientedCuboid<DIMENSION>> ();

	Shape<DIMENSION, 0>::setShapeStaticInfo(shapeInfo);
	Shape<DIMENSION, 0>::setEarlyRejectionEnabled(false);     // it doesn't help here
}

template <unsigned short DIMENSION>
bool OrientedCuboid<DIMENSION>::overlap(BoundaryConditions<DIMENSION> *bc, const Shape<DIMENSION, 0> *s) const{
    return this->pointInside(bc, s->getPosition(), s->getOrientation(), 0.0);
}

template <unsigned short DIMENSION>
bool OrientedCuboid<DIMENSION>::pointInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &position,
										   const Orientation<0> &orientation, double orientationRange) const{
	Vector<DIMENSION> ta = bc->getTranslation(this->getPosition(), position);
	for(unsigned short i=0; i<DIMENSION; i++){
		if (std::fabs(this->getPosition()[i] - (position[i] + ta[i])) > size[i]){
			return false;
		}
	}
	return true;
}



/*
template <unsigned short DIMENSION>
bool OrientedCuboid<DIMENSION>::voxelInside(BoundaryConditions<DIMENSION> *bc, Voxel *v, double spatialSize, double angularSize, const std::vector<const RSAShape *> *shapes){
	// check over single shapes
	RSAVector vpos = v->getPosition();
	for(const RSAShape* s : *shapes){
		if (s->voxelInside(bc, vpos, v->getOrientation(), spatialSize, angularSize) )
			return true;
	}
	// check over pairs
	char result = 0;
	bool fbreak = false;
	for (size_t i=0; i<shapes->size(); i++){
		const RSAShape* si = (*shapes)[i];
		RSAVector pi = si->getPosition();
		for (size_t j=i+1; j<shapes->size(); j++){
			const RSAShape* sj = (*shapes)[j];
			RSAVector pj = sj->getPosition();
			for (ushort k=0; k<DIMENSION; k++){
				if (std::fabs(pi[k] - pj[k]) > 2*size[k]){ // this pair does not exclude any voxel
					fbreak = true;
					break;
				}
				else if (std::fabs(pi[k] - pj[k]) > size[k]) // voxel may be at the boundary of exclusion zones
					// check if voxel is between shapes
					if (
							((std::fabs(vpos[k] - pi[k])<size[k] || std::fabs(vpos[k] - pj[k])<size[k])) &&
							((std::fabs(vpos[k] + spatialSize - pi[k])<size[k] || std::fabs(vpos[k] + spatialSize - pj[k])<size[k]))
					)
						result++;
					else{
						fbreak = true;
						break;
					}
				else  // std::fabs(pi[k] - pj[k]) < size[k]
					// another check
					if (!
							((std::fabs(vpos[k] - pi[k])<0.5*size[k] || std::fabs(vpos[k] - pj[k])<0.5*size[k])) &&
							((std::fabs(vpos[k] + spatialSize - pi[k])<0.5*size[k] || std::fabs(vpos[k] + spatialSize - pj[k])<0.5*size[k]))
					){
						fbreak = true;
						break;
					}
			}
			if (fbreak == false && result < DIMENSION)
				return true;
		}
	}
	return false;
}
*/

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
