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
OrientedCuboid<DIMENSION>::OrientedCuboid() : ConvexShape<DIMENSION, 0>(){
}

template <unsigned short DIMENSION>
OrientedCuboid<DIMENSION>::~OrientedCuboid() {
}

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
        
    // Calculate static params
//    volume = std::accumulate(size, size + staticDimension, 1.0, std::multiplies<double>());
 
   	Shape<DIMENSION, 0>::setNeighbourListCellSize(*std::max_element(OrientedCuboid<DIMENSION>::size, OrientedCuboid<DIMENSION>::size + DIMENSION));
   	Shape<DIMENSION, 0>::setVoxelSpatialSize(0.5*(*std::min_element(OrientedCuboid<DIMENSION>::size, OrientedCuboid<DIMENSION>::size + DIMENSION)));
//	Shape<DIMENSION, 0>::setCreateShapeImpl([](RND *rnd) -> Shape<DIMENSION, 0>* {
//		return new OrientedCuboid<DIMENSION>;
//	});
	Shape<DIMENSION, 0>::template setDefaultCreateShapeImpl <OrientedCuboid<DIMENSION>> ();
}

template <unsigned short DIMENSION>
bool OrientedCuboid<DIMENSION>::overlap(BoundaryConditions<DIMENSION> *bc, const Shape<DIMENSION, 0> *s) const{
    return this->pointInside(bc, s->getPosition());
}

template <unsigned short DIMENSION>
bool OrientedCuboid<DIMENSION>::pointInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &da) const{
//    if (OrientedCuboid<DIMENSION>::do2Drotation){
//    	return Cuboid::pointInside(bc, da);
//    }else{
    	Vector<DIMENSION> ta = bc->getTranslation(this->getPosition(), da);
    	for(unsigned short i=0; i<DIMENSION; i++){
    		if (std::fabs(this->getPosition()[i] - (da[i] + ta[i])) > size[i]){
    			return false;
    		}
    	}
    	return true;
//    }
}

template <unsigned short DIMENSION>
bool OrientedCuboid<DIMENSION>::pointInside(BoundaryConditions<DIMENSION> *bc, const Vector<DIMENSION> &position,
										   const Orientation<0> &orientation, double orientationRange) const{
	return 0;
}

template <unsigned short DIMENSION>
double OrientedCuboid<DIMENSION>::getVolume() const{
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