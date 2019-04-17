/*
 * Polygon.cpp
 *
 *  Created on: 16.04.2018
 *      Author: ciesla
 */

#include "Polygon.h"
#include <cmath>
#include <sstream>


std::vector<double> Polygon::vertexR;
std::vector<double> Polygon::vertexTheta;
std::vector<std::pair<size_t, size_t>> Polygon::segments;
std::vector<std::pair<size_t, size_t>> Polygon::helperSegments;
double Polygon::inscribedCircleRadius;

double Polygon::getCircumscribedCircleRadius(){
	double result = 0.0;
	for(std::pair<size_t, size_t> segment: Polygon::segments){
		if (result < Polygon::vertexR[segment.first])
			result = Polygon::vertexR[segment.first];
		if (result < Polygon::vertexR[segment.second])
			result = Polygon::vertexR[segment.second];
	}
	return result;
}

double Polygon::getInscribedCircleRadius(){
	if (Polygon::vertexR.size() < 3){
		return 0.0;
	}

	double result = Polygon::vertexR[Polygon::segments[0].first];
	for (size_t i = 0; i < Polygon::segments.size(); i++){
		std::pair<size_t, size_t> segment = Polygon::segments[i];
		double x3 = Polygon::vertexR[segment.first] * std::cos(Polygon::vertexTheta[segment.first]);
		double y3 = Polygon::vertexR[segment.first] * std::sin(Polygon::vertexTheta[segment.first]);
		double x4 = Polygon::vertexR[segment.second] * std::cos(Polygon::vertexTheta[segment.second]);
		double y4 = Polygon::vertexR[segment.second] * std::sin(Polygon::vertexTheta[segment.second]);
		double t = (x4 - x3) / (y4 - y3);
		double d = std::abs(t*y3 - x3) / std::sqrt(1.0 + t*t);
		result = std::min(result, d);
	}
	return result;
}

//calculate the area of the triangle made from the origin, vertex i, and vertex j
double Polygon::getTriangleArea(size_t i, size_t j){

	//a, b, and c are the side lengths of the triangle
	double a = Polygon::vertexR[i];
	double b = Polygon::vertexR[j];
	double dxc = Polygon::vertexR[i] * std::cos(Polygon::vertexTheta[i]) - Polygon::vertexR[j] * std::cos(Polygon::vertexTheta[j]);
	double dyc = Polygon::vertexR[i] * std::sin(Polygon::vertexTheta[i]) - Polygon::vertexR[j] * std::sin(Polygon::vertexTheta[j]);
	double c = std::sqrt(dxc*dxc + dyc*dyc);

	double s = 0.5*(a + b + c);
	return std::sqrt(s*(s - a)*(s - b)*(s - c));
}

/**
 * moves the center of the polygon to (0,0)
 */
void Polygon::centerPolygon(){
	double cx = 0.0, cy = 0.0;
	for (size_t i=0; i<Polygon::vertexR.size(); i++){
		cx += Polygon::vertexR[i]*std::cos(Polygon::vertexTheta[i]);
		cy += Polygon::vertexR[i]*std::sin(Polygon::vertexTheta[i]);
	}
	cx /= Polygon::vertexR.size();
	cy /= Polygon::vertexR.size();
	double x, y;
	for (size_t i=0; i<Polygon::vertexR.size(); i++){
		x = Polygon::vertexR[i]*std::cos(Polygon::vertexTheta[i]) - cx;
		y = Polygon::vertexR[i]*std::sin(Polygon::vertexTheta[i]) - cy;
		Polygon::vertexR[i] = std::sqrt(x*x + y*y);
		Polygon::vertexTheta[i] = std::atan2(y, x);
	}
}

void Polygon::createStarHelperSegments(){
	Polygon::helperSegments.clear();
	Polygon::vertexR.push_back(0.0);
	Polygon::vertexTheta.push_back(0.0);
	for (size_t i=0; i<Polygon::vertexR.size()-1; i++){
		std::pair<size_t, size_t> segment(i, Polygon::vertexR.size()-1);
		Polygon::helperSegments.push_back(segment);
	}
}

/**
 * @param args contains information about coordinates of each vertex of the polygon. Adjacent vertices are connected.
 * Data should be separated by only spaces, and should be: number_of_vertices [xy/rt] c01 c02 c11 c12 c21 c22 ...
 * xy means cartesian coordinates, and rt means polar coordinates
 * Example format of coordinates
 * 4 xy 1 1 1 -1 -1 -1 -1 1 4 0 1 1 2 2 3 3 0 0
 * or equivalently
 * 4 rt 1.4142 0.7854 1.4142 2.3562 1.4142 3.9270 1.4142 5.4978 4 0 1 1 2 2 3 3 0 0
 */
void Polygon::initClass(const std::string &args){
	std::istringstream in(args);

	size_t n;
	in >> n;
	std::string format;
	in >> format;
	double r, t;
	// reading vertices
	for (size_t i=0; i<n; i++){
		double c1, c2;
		in >> c1;
		in >> c2;
		if(format == "xy"){
			r = std::sqrt(c1*c1 + c2*c2);
			t = std::atan2(c2, c1);
		}else if (format == "rt"){
			r = c1;
			t = c2;
		}else{
			throw std::runtime_error("Wrong coordinate format. Use rt or xy");
		}
		Polygon::vertexR.push_back(r);
		Polygon::vertexTheta.push_back(t);
	}

	in >> n;
	size_t segBeg;
	in >> segBeg;
    if (segBeg > Polygon::vertexR.size()){
        throw std::runtime_error("Wrong vertex number in segment");
    }
	// reading segments
	for(size_t i = 1; i<n; i++){
		size_t segEnd;
		in >> segEnd;
		if (segEnd >Polygon::vertexR.size()){
			throw std::runtime_error("Wrong vertex number in segment");
		}
		Polygon::segments.emplace_back(segBeg, segEnd);
		segBeg = segEnd;
	}
	Polygon::segments.emplace_back(Polygon::segments.back().second, Polygon::segments.front().first);

	if (args.find("starHelperSegments")!=std::string::npos){
		Polygon::centerPolygon();
		Polygon::createStarHelperSegments();
	}else{
		in >> n;
		// reading helper segments
		for(size_t i = 0; i<n; i++){
			size_t i1, i2;
			in >> i1;
			in >> i2;
			if (i1>Polygon::vertexR.size() || i1>Polygon::vertexR.size()){
				throw std::runtime_error("Wrong vertex number in helper segment");
			}
			Polygon::helperSegments.emplace_back(i1, i2);
		}
	}

	double area = 0.0;
	for (auto segment : Polygon::segments)
	    area += Polygon::getTriangleArea(segment.first, segment.second);
	for (double & vR : Polygon::vertexR)
		vR /= std::sqrt(area);

	Shape<2, 1>::setNeighbourListCellSize(2.0*Polygon::getCircumscribedCircleRadius());
	Polygon::inscribedCircleRadius = Polygon::getInscribedCircleRadius();
	Shape<2, 1>::setVoxelSpatialSize(1.4*Polygon::inscribedCircleRadius);
	Shape<2, 1>::setVoxelAngularSize(2*M_PI);
	Shape<2, 1>::setSupportsSaturation(true);
	Shape<2, 1>::setDefaultCreateShapeImpl <Polygon> ();

	#ifdef CUDA_ENABLED
		Polygon::cuInit();
	#endif
}


	//test if line segment from point 1 to 2 intersects with line segment from point 3 to 4
bool Polygon::lineLineIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){
	double o1 = (y2 - y1)*(x3 - x2) - (x2 - x1)*(y3 - y2);
	double o2 = (y2 - y1)*(x4 - x2) - (x2 - x1)*(y4 - y2);
	double o3 = (y4 - y3)*(x1 - x4) - (x4 - x3)*(y1 - y4);
	double o4 = (y4 - y3)*(x2 - x4) - (x4 - x3)*(y2 - y4);

	return (o1*o2 < 0.0) && (o3*o4 < 0.0);
}

	//same as above, except that endpoints 3 and 4 comes from a line in a voxel, and thus carry an uncertainty
bool Polygon::lineVoxelIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double dx, double dtheta, double l3, double l4){
	double o1 = (y2 - y1)*(x3 - x2) - (x2 - x1)*(y3 - y2);
	double o2 = (y2 - y1)*(x4 - x2) - (x2 - x1)*(y4 - y2);
	double o3 = (y4 - y3)*(x1 - x4) - (x4 - x3)*(y1 - y4);
	double o4 = (y4 - y3)*(x2 - x4) - (x4 - x3)*(y2 - y4);

	if ((o1*o2 < 0.0) && (o3*o4 < 0.0)){
		//the voxel center is intersecting
		double d3 = dx + dtheta*l3;
		double d4 = dx + dtheta*l4;
		double do1 = d3*(std::abs(y2 - y1) + std::abs(x2 - x1));
		double do2 = d4*(std::abs(y2 - y1) + std::abs(x2 - x1));
		double do3 = (d4 + d3)*(std::abs(x1 - x4) + std::abs(y1 - y4) + 2 * d4) + d4*(std::abs(y4 - y3) + std::abs(x4 - x3));
		double do4 = (d4 + d3)*(std::abs(x2 - x4) + std::abs(y2 - y4) + 2 * d4) + d4*(std::abs(y4 - y3) + std::abs(x4 - x3));

		if (do1 < std::abs(o1) && do2 < std::abs(o2) && do3 < std::abs(o3) && do4 < std::abs(o4))
			return true;
		else
			return false;
	}
	else
		return false;
}

double Polygon::getVolume() const{
	double result = 0.0;
	for (size_t i = 0; i < Polygon::segments.size(); i++){
		std::pair<size_t, size_t> segment = Polygon::segments[i];
		result += Polygon::getTriangleArea(segment.first, segment.second);
	}
	return result;
}

#ifndef CUDA_ENABLED

bool Polygon::overlap(BoundaryConditions<2> *bc, const Shape<2, 1> *s) const{
	Polygon pol = dynamic_cast<const Polygon&>(*s);
	this->applyBC(bc, &pol);

	double polposition[2];
	pol.getPosition().copyToArray(polposition);
	double position[2];
	this->getPosition().copyToArray(position);

	//easy check
	double d2 = 0, tmp;
	for (unsigned short i = 0; i < 2; i++){
		tmp = position[i] - polposition[i];
		d2 += tmp*tmp;
	}
	if (std::sqrt(d2) < 2.0*Polygon::inscribedCircleRadius)
		return true;

	double angle = this->getOrientation()[0];
	double polangle = pol.getOrientation()[0];
	//complex check
	for (size_t i = 0; i < Polygon::segments.size() + Polygon::helperSegments.size(); i++){
		std::pair<size_t, size_t> polsegment;
		if (i<Polygon::segments.size()){
			polsegment = Polygon::segments[i];
		}else{
			polsegment = Polygon::helperSegments[i-Polygon::segments.size()];
		}
		double x1 = polposition[0] + Polygon::vertexR[polsegment.first] * std::cos(Polygon::vertexTheta[polsegment.first] + polangle);
		double y1 = polposition[1] + Polygon::vertexR[polsegment.first] * std::sin(Polygon::vertexTheta[polsegment.first] + polangle);
		double x2 = polposition[0] + Polygon::vertexR[polsegment.second] * std::cos(Polygon::vertexTheta[polsegment.second] + polangle);
		double y2 = polposition[1] + Polygon::vertexR[polsegment.second] * std::sin(Polygon::vertexTheta[polsegment.second] + polangle);

		for (size_t j = 0; j < Polygon::segments.size() + Polygon::helperSegments.size(); j++){
			std::pair<size_t, size_t> segment;
			if (j<Polygon::segments.size()){
				segment = Polygon::segments[j];
			}else{
				segment = Polygon::helperSegments[j-Polygon::segments.size()];
			}
			double x3 = position[0] + Polygon::vertexR[segment.first] * std::cos(Polygon::vertexTheta[segment.first] + angle);
			double y3 = position[1] + Polygon::vertexR[segment.first] * std::sin(Polygon::vertexTheta[segment.first] + angle);
			double x4 = position[0] + Polygon::vertexR[segment.second] * std::cos(Polygon::vertexTheta[segment.second] + angle);
			double y4 = position[1] + Polygon::vertexR[segment.second] * std::sin(Polygon::vertexTheta[segment.second] + angle);
			if (Polygon::lineLineIntersect(x1, y1, x2, y2, x3, y3, x4, y4))
				return true;
		}
	}
	return false;
}
#endif

bool Polygon::voxelInside(BoundaryConditions<2> *bc, const Vector<2> &voxelPosition,
						  const Orientation<1> &voxelOrientation, double spatialSize, double angularSize) const {
    if (voxelOrientation[0] > Shape<2, 1>::getVoxelAngularSize())
		return true;

	double halfSpatialSize = 0.5*spatialSize;
	Vector<2> voxelTranslation = bc->getTranslation(this->getPosition(), voxelPosition);
	Vector<2> spatialCenter = voxelPosition + voxelTranslation + Vector<2>{{halfSpatialSize, halfSpatialSize}};

    if (voxelInsideEasyCheck(spatialCenter, halfSpatialSize))
		return true;

    if (voxelInsideFullAngleCheck(spatialCenter, halfSpatialSize))
        return true;

    double halfAngularSize = 0.5*angularSize;
    double angularCenter = voxelOrientation[0] + halfAngularSize;

    return voxelInsideComplexCheck(spatialCenter, halfSpatialSize, angularCenter, halfAngularSize);
}

bool Polygon::voxelInsideEasyCheck(const Vector<2> &spatialCenter, double halfSpatialSize) const {
    Vector<2> voxelDistance = this->getPosition() - spatialCenter;
    for (unsigned short j = 0; j < 2; j++) {
        if (voxelDistance[j] > 0)
            voxelDistance[j] += halfSpatialSize;
        else
            voxelDistance[j] -= halfSpatialSize;
    }
    return voxelDistance.norm2() < 4*inscribedCircleRadius*inscribedCircleRadius;
}

bool Polygon::voxelInsideFullAngleCheck(const Vector<2> &spatialCenter, double halfSpatialSize) const {
    double pushDistance = Polygon::inscribedCircleRadius - halfSpatialSize * M_SQRT2;
    Assert(pushDistance >= 0);

    return this->pointInsidePushed(spatialCenter, pushDistance);
}

bool Polygon::pointInsidePushed(const Vector<2> &point, double pushDistance) const {
    if (this->pointInsidePolygon(point))
        return true;

    return false;
}

bool Polygon::pointInsidePolygon(const Vector<2> &point) const {
    std::vector<double> segmentXIntersections;

    for (auto segment : Polygon::segments) {
        Vector<2> s1 = Polygon::getVertexPosition(segment.first);
        Vector<2> s2 = Polygon::getVertexPosition(segment.second);

        if (s1[1] > s2[1])
            std::swap(s1, s2);

        if (point[1] > s2[1] || point[1] < s1[1])
            continue;

        double xIntersection = (s2[0] * (point[1] - s1[1]) - s2[0] * (point[1] - s1[1])) / (s2[1] - s1[1]);
        segmentXIntersections.push_back(xIntersection);
    }

    double px = point[0];
    auto onLeft = [px](double x0) { return x0 < px; };
    std::size_t intersectionsOnLeft = std::count_if(segmentXIntersections.begin(), segmentXIntersections.end(), onLeft);

    return intersectionsOnLeft % 2 == 1;
}

bool Polygon::voxelInsideComplexCheck(const Vector<2> &spatialCenter, double halfSpatialSize, double angularCenter,
                                      double halfAngularSize) const {
    Vector<2> position = this->getPosition();
    double angle = getOrientation()[0];
    for (size_t i = 0; i < segments.size() + helperSegments.size(); i++) {
        std::pair<size_t, size_t> segment;
        if (i < segments.size()) {
            segment = segments[i];
        } else {
            segment = helperSegments[i - segments.size()];
        }

        double x1 = position[0] + vertexR[segment.first] * cos(vertexTheta[segment.first] + angle);
        double y1 = position[1] + vertexR[segment.first] * sin(vertexTheta[segment.first] + angle);
        double x2 = position[0] + vertexR[segment.second] * cos(vertexTheta[segment.second] + angle);
        double y2 = position[1] + vertexR[segment.second] * sin(vertexTheta[segment.second] + angle);

        for (size_t j = 0; j < segments.size() + helperSegments.size(); j++) {
            std::pair<size_t, size_t> vsegment;
            if (j < segments.size()) {
                vsegment = segments[j];
            } else {
                vsegment = helperSegments[j - segments.size()];
            }
            double x3 = spatialCenter[0] + vertexR[vsegment.first] * cos(vertexTheta[vsegment.first] + angularCenter);
            double y3 = spatialCenter[1] + vertexR[vsegment.first] * sin(vertexTheta[vsegment.first] + angularCenter);
            double x4 = spatialCenter[0] + vertexR[vsegment.second] * cos(vertexTheta[vsegment.second] + angularCenter);
            double y4 = spatialCenter[1] + vertexR[vsegment.second] * sin(vertexTheta[vsegment.second] + angularCenter);
            if (lineVoxelIntersect(x1, y1, x2, y2, x3, y3, x4, y4, halfSpatialSize, halfAngularSize,
                                   vertexR[vsegment.first], vertexR[vsegment.second])) {
                return true;
            }
        }
    }
    return false;
}

Shape<2, 1> *Polygon::clone() const {
    return new Polygon(*this);
}

std::string Polygon::toPovray() const{
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);
	out << "  polygon {" << Polygon::segments.size()+1 << ", ";
	this->vertexToPovray(Polygon::segments[0].first, out);
	for (size_t i=0; i < Polygon::segments.size(); i++) {
        this->vertexToPovray(Polygon::segments[i].second, out);
        if (i<Polygon::segments.size()-1)
        	out << " ,";
    }
	out << "  texture { finish { ambient 1 diffuse 0 } pigment { color Red} } }" << std::endl;

	return out.str();
}

std::string Polygon::toString() const {
	return this->toWolfram();
}

std::string Polygon::toWolfram() const {
    std::ostringstream out;
    out.precision(std::numeric_limits<double>::max_digits10);
    out << "Polygon[{";
    out << this->getVertexPosition(Polygon::segments[0].first);
    for (std::size_t i = 1; i < segments.size(); i++)
        out << ", " << this->getVertexPosition(Polygon::segments[i].first);
    out << "}]";
    return out.str();
}

Vector<2> Polygon::getVertexPosition(std::size_t index) const {
    Vector<2> position = this->getPosition();
    double angle = Polygon::vertexTheta[index] + this->getOrientation()[0];
    return Vector<2>{{
        position[0] + Polygon::vertexR[index] * std::cos(angle),
        position[1] + Polygon::vertexR[index] * std::sin(angle)
    }};
}

void Polygon::vertexToPovray(std::size_t index, std::ostream &out) const {
    Vector<2> position = this->getVertexPosition(index);
    out << "< " << position[0] << ", " << position[1] << ", 0.0002>";
}
