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
double OrientedCuboid<DIMENSION>::neighbourListCellSize;

template <unsigned short DIMENSION>
double OrientedCuboid<DIMENSION>::voxelSize;


template <unsigned short DIMENSION>
OrientedCuboid<DIMENSION>::OrientedCuboid() : Shape<DIMENSION>(){
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
 
   	neighbourListCellSize = *std::max_element(OrientedCuboid<DIMENSION>::size, OrientedCuboid<DIMENSION>::size + DIMENSION);
   	voxelSize = 0.5*(*std::min_element(OrientedCuboid<DIMENSION>::size, OrientedCuboid<DIMENSION>::size + DIMENSION));
}


template <unsigned short DIMENSION>
Shape<DIMENSION> * OrientedCuboid<DIMENSION>::create(RND *rnd){
	OrientedCuboid<DIMENSION> *cuboid;
//	if (OrientedCuboid::do2Drotation)
//	    cuboid = new OrientedCuboid(Matrix<3, 3>::rotation(
//	        0,
//	        std::asin(rnd->nextValue() * 2 - 1),
//	        rnd->nextValue() * 2 * M_PI));
//	else
    	cuboid = new OrientedCuboid<DIMENSION>();

#ifdef CUBOID_DEBUG
    std::cout << "Creating OrientedCuboid:" << std::endl;
    std::cout << cuboid->orientation;
#endif

    return cuboid;
}

template <unsigned short DIMENSION>
int OrientedCuboid<DIMENSION>::overlap(BoundaryConditions *bc, Shape<DIMENSION> *s){
    return this->pointInside(bc, s->getPosition());
}

template <unsigned short DIMENSION>
int OrientedCuboid<DIMENSION>::pointInside(BoundaryConditions *bc, double* da){
//    if (OrientedCuboid<DIMENSION>::do2Drotation){
//    	return Cuboid::pointInside(bc, da);
//    }else{
    	double ta[DIMENSION];
   		bc->getTranslation(ta, this->position, da);
    	for(unsigned short i=0; i<DIMENSION; i++){
    		if (std::fabs(this->position[i] - (da[i] + ta[i])) > size[i]){
    			return false;
    		}
    	}
    	return true;
//    }
}

template <unsigned short DIMENSION>
double OrientedCuboid<DIMENSION>::getNeighbourListCellSize(){
	return OrientedCuboid<DIMENSION>::neighbourListCellSize;
}

template <unsigned short DIMENSION>
double OrientedCuboid<DIMENSION>::getVoxelSize(){
	return OrientedCuboid<DIMENSION>::voxelSize;
}

template <unsigned short DIMENSION>
double OrientedCuboid<DIMENSION>::getVolume(){
	return std::accumulate(size, size + DIMENSION, 1.0, std::multiplies<double>());
}

template <unsigned short DIMENSION>
std::string OrientedCuboid<DIMENSION>::toPovray() const
{	
	std::string s = "";
	if (DIMENSION==2){
		double factor = 0.5;
		s += "  polygon {4, ";

		s += "< " + std::to_string(this->position[0] - factor*size[0]) + ", ";
		s += 		std::to_string(this->position[1] - factor*size[1]) + ", 0.1>, ";

		s += "< " + std::to_string(this->position[0] - factor*size[0]) + ", ";
		s += 		std::to_string(this->position[1] + factor*size[1]) + ", 0.1>, ";

		s += "< " + std::to_string(this->position[0] + factor*size[0]) + ", ";
		s += 		std::to_string(this->position[1] + factor*size[1]) + ", 0.1>, ";

		s += "< " + std::to_string(this->position[0] + factor*size[0]) + ", ";
		s += 		std::to_string(this->position[1] - factor*size[1]) + ", 0.1> ";

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
			s += std::to_string(this->position[i]);
			if (i<DIMENSION-1)
				s+= ", ";
		}
		s += "\n    >\n";
		s += "    texture { pigment { color Red } }\n  }\n";
	}
	return s;
}