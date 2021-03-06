/*
 * Polydisk.cpp
 *
 *  Created on: 23.09.2018
 *      Author: ciesla
 */

#include "Polydisk.h"
#include "../../../RND.h"
#include <cmath>
#include <sstream>


std::vector<double> Polydisk::diskCentreR;
std::vector<double> Polydisk::diskCentreTheta;
std::vector<double> Polydisk::diskR;

/* If one wants to center automatically on the largest disk:
 * number_of_disks [xy|rt] c01 c02 r0 c11 c12 r1 c21 c22 r2 ... area
 *
 * If one wants to retain original centering:
 * number_of_disks [xy|rt] c01 c02 r0 c11 c12 r1 c21 c22 r2 ... area 'dontcenter' insphere_radius */
void Polydisk::initClass(const std::string &args){
    Polydisk::diskCentreR.clear();
    Polydisk::diskCentreTheta.clear();
    Polydisk::diskR.clear();

	std::istringstream in(args);

	size_t n;
	in >> n;
	std::string format;
	in >> format;
	double r, t;
	for (size_t i=0; i<n; i++){
		double c1, c2, diskR_;
		in >> c1;
		in >> c2;
		in >> diskR_;
		if (diskR_ <= 0.0)
		    throw std::runtime_error("diskR <= 0 for " + std::to_string(i) + " coord");

		if(format == "xy"){
			r = std::sqrt(c1*c1 + c2*c2);
			t = std::atan2(c2, c1);
		}else if (format == "rt"){
			r = c1;
			t = c2;
		}else{
			throw std::runtime_error("Wrong coordinate format. Use rt or xy");
		}
		Polydisk::diskCentreR.push_back(r);
		Polydisk::diskCentreTheta.push_back(t);
		Polydisk::diskR.push_back(diskR_);
	}

	double initialArea;
	in >> initialArea;

	if (initialArea==0.0){
        initialArea = Polydisk::mcArea(100000000);
	}
	Polydisk::normalizeArea(initialArea);

	std::string dontCenterPhrase;
	in >> dontCenterPhrase;
	double insphereRadius{};
	if (dontCenterPhrase == "dontcenter") {
	    in >> insphereRadius;
	    Validate(insphereRadius > 0);
	    insphereRadius /= std::sqrt(initialArea);
	} else {
        Vector<2> newCentre = Polydisk::findPositionOfLargestDiskClosestToMassCentre();
        Polydisk::centerPolydisk(newCentre);
        insphereRadius = Polydisk::diskR[Polydisk::getIndicesOfLargestDisks().front()];
	}

    double circumsphereRadius = 0;
    for (size_t i = 0; i < Polydisk::diskCentreR.size(); i++)
        if (circumsphereRadius < Polydisk::diskCentreR[i] + Polydisk::diskR[i])
            circumsphereRadius = Polydisk::diskCentreR[i] + Polydisk::diskR[i];
    Validate(insphereRadius <= circumsphereRadius);

	ShapeStaticInfo<2, 1> shapeInfo;

	shapeInfo.setCircumsphereRadius(circumsphereRadius);
	shapeInfo.setInsphereRadius(insphereRadius);
	shapeInfo.setAngularVoxelSize(2*M_PI);
	shapeInfo.setSupportsSaturation(true);
	shapeInfo.setDefaultCreateShapeImpl <Polydisk> ();

	Shape::setShapeStaticInfo(shapeInfo);
}

void Polydisk::centerPolydisk(const Vector<2> &newCentre) {
    for (std::size_t i = 0; i < diskCentreR.size(); i++) {
        Vector<2> position = Polydisk::getStaticDiskPosition(i);
        position -= newCentre;

        Polydisk::diskCentreR[i] = position.norm();
        Polydisk::diskCentreTheta[i] = position.angle();
    }
}

Vector<2> Polydisk::findPositionOfLargestDiskClosestToMassCentre() {
    std::vector<std::size_t> indicesOfLargestDisks = Polydisk::getIndicesOfLargestDisks();

    std::vector<Vector<2>> positionsOfLargestDisks;
    std::transform(indicesOfLargestDisks.begin(), indicesOfLargestDisks.end(),
                   std::back_inserter(positionsOfLargestDisks), &Polydisk::getStaticDiskPosition);

    Vector<2> approximateMassCentre = Polydisk::getApproximateMassCentre();

    return *std::min_element(positionsOfLargestDisks.begin(), positionsOfLargestDisks.end(),
                             [approximateMassCentre](auto p1, auto p2) {
                                 return (p1 - approximateMassCentre).norm2() < (p2 - approximateMassCentre).norm2();
                             });
}

Vector<2> Polydisk::getApproximateMassCentre() {
    Vector<2> approximateMassCentre;
    for (size_t i = 0; i < diskCentreR.size(); i++)
        approximateMassCentre += getStaticDiskPosition(i) * std::pow(diskR[i], 2);
    double sumOfDiskRadiusSquares = std::accumulate(diskR.begin(), diskR.end(), 0.,
                                                    [](double sum, double r) { return sum + r*r; });
    approximateMassCentre /= sumOfDiskRadiusSquares;
    return approximateMassCentre;
}

std::vector<std::size_t> Polydisk::getIndicesOfLargestDisks() {
    std::vector<std::size_t> indexesOfLargestDisks = {0};
    double radiusOfLargestDisk = diskR[0];
    for (std::size_t i = 1; i < diskR.size(); i++) {
        double diskR = Polydisk::diskR[i];
        if (diskR == radiusOfLargestDisk) {
            indexesOfLargestDisks.push_back(i);
        } else if (diskR > radiusOfLargestDisk) {
            radiusOfLargestDisk = diskR;
            indexesOfLargestDisks = {i};
        }
    }
    return indexesOfLargestDisks;
}

double Polydisk::mcArea(size_t mcTrials){
	double maxSize = 1.0;
	for (size_t i = 0; i < Polydisk::diskCentreR.size(); i++){
		if (maxSize < Polydisk::diskCentreR[i] + Polydisk::diskR[i])
			maxSize = Polydisk::diskCentreR[i] + Polydisk::diskR[i];
	}
	RND rnd(12345);
	double x, y, dx, dy;
	size_t inside = 0;

    for (size_t i=0; i<mcTrials; i++){
		x = (2.0*rnd.nextValue()-1.0)*maxSize;
		y = (2.0*rnd.nextValue()-1.0)*maxSize;
		for(size_t disk = 0; disk < Polydisk::diskR.size(); disk++){
			dx = x - Polydisk::diskCentreR[disk]*std::cos(Polydisk::diskCentreTheta[disk]);
			dy = y - Polydisk::diskCentreR[disk]*std::sin(Polydisk::diskCentreTheta[disk]);
			if (dx*dx + dy*dy < Polydisk::diskR[disk]*Polydisk::diskR[disk]){
				inside++;
				break;
			}
		}
	}
	return ((double)inside/(double)mcTrials)*4.0*maxSize*maxSize;
}

void Polydisk::normalizeArea(double area){
	double denominator = std::sqrt(area);

	for(size_t disk = 0; disk < Polydisk::diskR.size(); disk++){
		Polydisk::diskCentreR[disk] /= denominator;
		Polydisk::diskR[disk] /= denominator;
	}
}

bool Polydisk::diskDiskIntersect(size_t disk0, const Vector<2> shape0Position, const Orientation<1> shape0Orientation,
								 size_t disk1, const Vector<2> shape1Position, const Orientation<1> shape1Orientation){
	double x0 = shape0Position[0] + Polydisk::diskCentreR[disk0]*std::cos(Polydisk::diskCentreTheta[disk0] + shape0Orientation[0]);
	double y0 = shape0Position[1] + Polydisk::diskCentreR[disk0]*std::sin(Polydisk::diskCentreTheta[disk0] + shape0Orientation[0]);

	double x1 = shape1Position[0] + Polydisk::diskCentreR[disk1]*std::cos(Polydisk::diskCentreTheta[disk1] + shape1Orientation[0]);
	double y1 = shape1Position[1] + Polydisk::diskCentreR[disk1]*std::sin(Polydisk::diskCentreTheta[disk1] + shape1Orientation[0]);

	double distance2 = (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1);
	distance2 -= (Polydisk::diskR[disk0] + Polydisk::diskR[disk1])*(Polydisk::diskR[disk0] + Polydisk::diskR[disk1]);
	return distance2 < 0.0;
}

Vector<2> Polydisk::minmaxCos(double theta, double dt){

	while (theta>2*M_PI)
		theta -= 2*M_PI;
	while (theta<0.0)
		theta += 2*M_PI;

	double v1 = std::cos(theta);
	double v2 = std::cos(theta+dt);
	Vector<2> result{{
		std::min(v1, v2),
		std::max(v1,v2)
    }};

	if ( (theta < 0 && theta+dt > 0 ) || (theta < 2.0*M_PI && theta+dt > 2.0*M_PI) ){
		result[1] = 1.0;
	}else if (theta < M_PI && theta+dt > M_PI){
		result[0] = -1.0;
	}
	return result;
}

Vector<2> Polydisk::minmaxSin(double theta, double dt){

	while (theta>2*M_PI)
		theta -= 2*M_PI;
	while (theta<0.0)
		theta += 2*M_PI;

	double v1 = std::sin(theta);
	double v2 = std::sin(theta+dt);

	Vector<2> result{{
		std::min(v1, v2),
		std::max(v1,v2)
    }};
	if (theta < 0.5*M_PI && theta+dt > 0.5*M_PI ){
		result[1] = 1.0;
	}else if (theta < 1.5*M_PI && theta+dt > 1.5*M_PI){
		result[0] = -1.0;
	}
	return result;
}

// checks if disk1 beeing a part of a particle inside the voxel intersects with disk0 beeing a part of a particle in packing
bool Polydisk::diskVoxelIntersect(size_t disk0, const Vector<2> shape0Position, const Orientation<1> shape0Orientation,
		 	 	 	 	 	 	  size_t disk1, const Vector<2> shape1Position, const Orientation<1> shape1Orientation,
								  double spatialSize, double angularSize){

	Vector<2> minmax = Polydisk::minmaxCos(Polydisk::diskCentreTheta[disk1]+shape1Orientation[0], angularSize);
	// minmaxX stores ranges for x coordinate on the disk1 (inside voxel)
	Vector<2> minmaxX = Vector<2>{{
		shape1Position[0] + Polydisk::diskCentreR[disk1]*minmax[0],
		shape1Position[0] + Polydisk::diskCentreR[disk1]*minmax[1] + spatialSize
	}};
	minmax = Polydisk::minmaxSin(Polydisk::diskCentreTheta[disk1]+shape1Orientation[0], angularSize);
	Vector<2> minmaxY = Vector<2>{{
		shape1Position[1] + Polydisk::diskCentreR[disk1]*minmax[0],
		shape1Position[1] + Polydisk::diskCentreR[disk1]*minmax[1] + spatialSize
	}};

	double maxDistance2 = 0.0;
	double x0 = shape0Position[0] + Polydisk::diskCentreR[disk0]*std::cos(Polydisk::diskCentreTheta[disk0]+shape0Orientation[0]);
	if (x0 < minmaxX[0])
		maxDistance2 += (x0-minmaxX[1])*(x0-minmaxX[1]);
	else if (x0 >= minmaxX[0] && x0 <= minmaxX[1])
		maxDistance2 += std::max((x0-minmaxX[0])*(x0-minmaxX[0]), (x0-minmaxX[1])*(x0-minmaxX[1]));
	else // disk0[0] > x1minmax[0]
		maxDistance2 += (x0-minmaxX[0])*(x0-minmaxX[0]);

	double y0 = shape0Position[1] + Polydisk::diskCentreR[disk0]*std::sin(Polydisk::diskCentreTheta[disk0]+shape0Orientation[0]);
	if (y0 < minmaxY[0])
		maxDistance2 += (y0-minmaxY[1])*(y0-minmaxY[1]);
	else if (y0 >= minmaxY[0] && y0 <= minmaxY[1])
		maxDistance2 += std::max((y0-minmaxY[0])*(y0-minmaxY[0]), (y0-minmaxY[1])*(y0-minmaxY[1]));
	else // disk0[0] > x1minmax[0]
		maxDistance2 += (y0-minmaxY[0])*(y0-minmaxY[0]);

	maxDistance2 -= (Polydisk::diskR[disk0] + Polydisk::diskR[disk1])*(Polydisk::diskR[disk0] + Polydisk::diskR[disk1]);
	return (maxDistance2<0.0);
}


double Polydisk::getVolume(unsigned short dim) const {
    if (dim != 2)
        throw std::runtime_error ("Polydisk supports only 2D packings");

    return 1.0;
}

bool Polydisk::overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const{
    switch (this->overlapEarlyRejection(bc, s)) {
        case TRUE:      return true;
        case FALSE:     return false;
        case UNKNOWN:   break;
    }

	Polydisk pol = dynamic_cast<const Polydisk&>(*s);
	this->applyBC(bc, &pol);

	Vector<2> polPosition = pol.getPosition();
	Orientation<1> polOrientation = pol.getOrientation();

	Vector<2> thisPosition = this->getPosition();
	Orientation<1> thisOrientation = this->getOrientation();

	for (size_t i = 0; i < Polydisk::diskCentreR.size(); i++){
		for (size_t j = 0; j < Polydisk::diskCentreR.size(); j++){
			if ( Polydisk::diskDiskIntersect(i, polPosition, polOrientation, j, thisPosition, thisOrientation) )
				return true;
		}
	}
	return false;
}

bool Polydisk::fullAngleVoxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition,
                                    double spatialSize) const
{
    double voxelCircumsphereRadius = spatialSize / M_SQRT2;
    double insphereRadius = Polydisk::getInsphereRadius();

    for (std::size_t i = 0; i < Polydisk::diskR.size(); i++) {
        double voxelDistance2 = bc->distance2(this->getDiskPosition(i), voxelPosition);
        if (voxelDistance2 <= std::pow(Polydisk::diskR[i] + insphereRadius - voxelCircumsphereRadius, 2))
            return true;
    }
    return false;
}

bool Polydisk::voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition, const Orientation<1> &voxelOrientation, double spatialSize, double angularSize) const{

	if (voxelOrientation[0] > Shape<2, 1>::getAngularVoxelSize())
		return true;

    switch(this->voxelInsideEarlyRejection(bc, voxelPosition, voxelOrientation, spatialSize, angularSize)) {
        case TRUE:      return true;
        case FALSE:     return false;
        case UNKNOWN:   break;
    }

    // Version optimized for full angle - it was made mainly for OrientedFibrinogen, because it took ages to generate
    if (angularSize >= 2*M_PI)
	    return this->fullAngleVoxelInside(bc, voxelPosition, spatialSize);

	Vector<2> translation = bc->getTranslation(voxelPosition, this->getPosition());
	Vector<2> thisPosition = this->getPosition() + translation;
	Orientation<1> thisOrientation = this->getOrientation();

	// loop over disks in this particle
	for (size_t i = 0; i < Polydisk::diskCentreR.size(); i++){
		// loop over disks in virtual particle inside voxel
		for (size_t j = 0; j < Polydisk::diskCentreR.size(); j++){
			// if a disk of virtual particle (placed anywhere in the voxel) intersects with disk of this particle, the voxel is fully inside the exclusion zone
			if (Polydisk::diskVoxelIntersect(i, thisPosition, thisOrientation, j, voxelPosition, voxelOrientation, spatialSize, angularSize) )
				return true;
		}
	}
	return false;
}

Shape<2, 1> *Polydisk::clone() const {
    return new Polydisk(*this);
}

std::string Polydisk::toPovray() const{
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);
	for (size_t i=0; i < Polydisk::diskCentreR.size(); i++) {
        Vector<2> diskPosition = this->getDiskPosition(i);
		double r = Polydisk::diskR[i];
		out << "  disc { < " << diskPosition[0] << ", " << diskPosition[1] << ", 0.05>, <0.0, 0.0, 1.0>, ";
		out << r << std::endl;
		out << "    texture { finish { ambient 1 diffuse 0 } pigment { color Red} }" << std::endl << "}" << std::endl;
    }

	return out.str();
}

std::string Polydisk::toWolfram() const {
    std::stringstream out;
    out.precision(std::numeric_limits<double>::max_digits10);

    out << std::fixed << "{";
    for (std::size_t i = 0, max = Polydisk::diskCentreR.size(); i < max; i++) {
        out << "Disk[" << this->getDiskPosition(i) << ", " << Polydisk::diskR[i] << "]";
        if (i < max - 1)
            out << ", ";
    }
    out << "}";

    return out.str();
}

Vector<2> Polydisk::getDiskPosition(std::size_t diskIndex) const {
    Expects(diskIndex < Polydisk::diskCentreR.size());

    auto position = this->getPosition();
    double angle = this->getOrientation()[0];
    return {{position[0] + Polydisk::diskCentreR[diskIndex] * std::cos(Polydisk::diskCentreTheta[diskIndex] + angle),
                    position[1] + Polydisk::diskCentreR[diskIndex] * std::sin(Polydisk::diskCentreTheta[diskIndex] + angle)}};
}


Vector<2> Polydisk::getStaticDiskPosition(std::size_t diskIndex) {
    Expects(diskIndex < Polydisk::diskCentreR.size());

    return {{Polydisk::diskCentreR[diskIndex] * std::cos(Polydisk::diskCentreTheta[diskIndex]),
             Polydisk::diskCentreR[diskIndex] * std::sin(Polydisk::diskCentreTheta[diskIndex])}};
}

