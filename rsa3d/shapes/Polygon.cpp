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
double Polygon::inscribedCircleRadius;

double Polygon::getCircumscribedCircleRadius(){
	double result = 0.0;
	for(double r: Polygon::vertexR){
		if (result < r)
			result = r;
	}
	return result;
}

double Polygon::getInscribedCircleRadius(){
	if (Polygon::vertexR.size() < 3){
		return 0.0;
	}

	double result = Polygon::vertexR[0];
	for (size_t i = 0; i < Polygon::vertexR.size(); i++){
		size_t next = (i + 1) % Polygon::vertexR.size();
		double x3 = Polygon::vertexR[i] * std::cos(Polygon::vertexTheta[i]);
		double y3 = Polygon::vertexR[i] * std::sin(Polygon::vertexTheta[i]);
		double x4 = Polygon::vertexR[next] * std::cos(Polygon::vertexTheta[next]);
		double y4 = Polygon::vertexR[next] * std::sin(Polygon::vertexTheta[next]);
		double t = (x4 - x3) / (y4 - y3);
		double d = std::abs(t*y3 - x3) / std::sqrt(1.0 + t*t);
		result = std::min(result, d);
	}
	return result;
}

//calculate the area of the triangle made from the origin, vertex i, and vertex (i+1)
double Polygon::getTriangleArea(size_t i){
	size_t j = (i + 1) % Polygon::vertexR.size();

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
 * @param args contains information about coordinates of each vertex of the polygon. Adjacent vertices are connected.
 * Data should be separated by only spaces, and should be: number_of_vertices [xy/rt] c01 c02 c11 c12 c21 c22 ...
 * xy means cartesian coordinates, and rt means polar coordinates
 * Example format of coordinates
 * 4 xy 1 1 1 -1 -1 -1 -1 1
 * or equivalently
 * 4 rt 1.4142 0.7854 1.4142 2.3562 1.4142 3.9270 1.4142 5.4978
 */
void Polygon::initClass(const std::string &args){
	std::istringstream in(args);

	size_t n;
	in >> n;
	std::string format;
	in >> format;
	double r, t;
	for (size_t i=0; i<n; i++){
		double c1, c2;
		in >> c1;
		in >> c2;
		if(format.compare("xy")==0){
			r = std::sqrt(c1*c1 + c2*c2);
			t = std::atan2(c2, c1);
		}else if (format.compare("rt")==0){
			r = c1;
			t = c2;
		}else{
			throw std::runtime_error("Wrong coordinate format. Use rt or xy");
		}
		Polygon::vertexR.push_back(r);
		Polygon::vertexTheta.push_back(t);
	}

	double area = 0.0;
	for (size_t i = 0; i < Polygon::vertexR.size(); i++){
		area += Polygon::getTriangleArea(i);
	}

	for(size_t i = 0; i<Polygon::vertexR.size(); i++){
		Polygon::vertexR[i] /= std::sqrt(area);
	}

	Shape<2, 1>::setNeighbourListCellSize(2.0*Polygon::getCircumscribedCircleRadius());
	Shape<2, 1>::setNeighbourListCellSize(2*Polygon::getCircumscribedCircleRadius());
	Polygon::inscribedCircleRadius = Polygon::getInscribedCircleRadius();
	Shape<2, 1>::setVoxelSpatialSize(1.4*Polygon::inscribedCircleRadius);
	Shape<2, 1>::setVoxelAngularSize(2*M_PI);
	Shape<2, 1>::setDefaultCreateShapeImpl <Polygon> ();
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

Polygon::Polygon() : AnisotropicShape2D(){
	this->setAngle(0.0);
}

double Polygon::getVolume(){
	double result = 0.0;
	for (size_t i = 0; i < Polygon::vertexR.size(); i++){
		result += Polygon::getTriangleArea(i);
	}
	return result;
}

bool Polygon::overlap(BoundaryConditions *bc, Shape<2, 1> *s) const{
	Polygon *polPtr = (Polygon *)s;
	Polygon pol(*polPtr);
	this->applyBC(bc, &pol);

	double *polposition = pol.getPosition();
	double *position = this->getPosition();

	//easy check
	double d2 = 0, tmp;
	for (unsigned short i = 0; i < RSA_SPATIAL_DIMENSION; i++){
		tmp = position[i] - polposition[i];
		d2 += tmp*tmp;
	}
	if (std::sqrt(d2) < 2.0*Polygon::inscribedCircleRadius)
		return true;

	double angle = this->getAngle();
	double polangle = pol.getAngle();
	//complex check
	for (size_t i = 0; i<Polygon::vertexR.size(); i++){
		size_t nexti = (i + 1) % Polygon::vertexR.size();
		double x1 = polposition[0] + Polygon::vertexR[i] * std::cos(Polygon::vertexTheta[i] + polangle);
		double y1 = polposition[1] + Polygon::vertexR[i] * std::sin(Polygon::vertexTheta[i] + polangle);
		double x2 = polposition[0] + Polygon::vertexR[nexti] * std::cos(Polygon::vertexTheta[nexti] + polangle);
		double y2 = polposition[1] + Polygon::vertexR[nexti] * std::sin(Polygon::vertexTheta[nexti] + polangle);

		for (size_t j = 0; j < Polygon::vertexR.size(); j++){
			size_t nextj = (j + 1) % Polygon::vertexR.size();
			double x3 = position[0] + Polygon::vertexR[j] * std::cos(Polygon::vertexTheta[j] + angle);
			double y3 = position[1] + Polygon::vertexR[j] * std::sin(Polygon::vertexTheta[j] + angle);
			double x4 = position[0] + Polygon::vertexR[nextj] * std::cos(Polygon::vertexTheta[nextj] + angle);
			double y4 = position[1] + Polygon::vertexR[nextj] * std::sin(Polygon::vertexTheta[nextj] + angle);
			if (Polygon::lineLineIntersect(x1, y1, x2, y2, x3, y3, x4, y4))
				return true;
		}
	}
	return false;
}

bool Polygon::voxelInside(BoundaryConditions *bc, const double *voxelPosition, const std::array<double, RSA_ANGULAR_DIMENSION> &voxelOrientation, double spatialSize, double angularSize) const{


	if (voxelOrientation[0] > Shape<2, 1>::getVoxelAngularSize())
		return true;

	double *position = this->getPosition();
	double translation[RSA_SPATIAL_DIMENSION];
	bc->getTranslation(translation, voxelPosition, position);

	double spatialCenter[RSA_SPATIAL_DIMENSION];
	double halfSpatialSize = 0.5*spatialSize;
	for(unsigned short i = 0; i<RSA_SPATIAL_DIMENSION; i++){
		spatialCenter[i] = voxelPosition[i] + halfSpatialSize;
	}

	double halfAngularSize = 0.5*angularSize;
	double angularCenter = voxelOrientation[0] + halfAngularSize;

	double angle = this->getAngle();
	//complex check
	for (size_t i = 0; i < Polygon::vertexR.size(); i++){
		size_t nexti = (i + 1) % Polygon::vertexR.size();
		double x1 = position[0] + translation[0] + Polygon::vertexR[i] * std::cos(Polygon::vertexTheta[i] + angle);
		double y1 = position[1] + translation[1] + Polygon::vertexR[i] * std::sin(Polygon::vertexTheta[i] + angle);
		double x2 = position[0] + translation[0] + Polygon::vertexR[nexti] * std::cos(Polygon::vertexTheta[nexti] + angle);
		double y2 = position[1] + translation[1] + Polygon::vertexR[nexti] * std::sin(Polygon::vertexTheta[nexti] + angle);

		for (size_t j = 0; j < Polygon::vertexR.size(); j++){
			size_t nextj = (j + 1) % Polygon::vertexR.size();
			double x3 = spatialCenter[0] + Polygon::vertexR[j] * std::cos(Polygon::vertexTheta[j] + angularCenter);
			double y3 = spatialCenter[1] + Polygon::vertexR[j] * std::sin(Polygon::vertexTheta[j] + angularCenter);
			double x4 = spatialCenter[0] + Polygon::vertexR[nextj] * std::cos(Polygon::vertexTheta[nextj] + angularCenter);
			double y4 = spatialCenter[1] + Polygon::vertexR[nextj] * std::sin(Polygon::vertexTheta[nextj] + angularCenter);
			if (Polygon::lineVoxelIntersect(x1, y1, x2, y2, x3, y3, x4, y4, halfSpatialSize, halfAngularSize, Polygon::vertexR[j], Polygon::vertexR[nextj]))
				return true;
		}
	}
	return false;
}

Shape<2, 1> *Polygon::clone() const {
    return new Polygon(*this);
}

std::string Polygon::toPovray() const{
	std::stringstream out;
	double *position = this->getPosition();
	out.precision(std::numeric_limits< double >::max_digits10);
	out << "  polygon {" << Polygon::vertexR.size()+1 << ", ";
	for (size_t i=0; i < Polygon::vertexR.size(); i++){
		out << "< ";
		out << position[0] + Polygon::vertexR[i] * std::cos(Polygon::vertexTheta[i] + this->getAngle());
		out << ", ";
		out << position[1] + Polygon::vertexR[i] * std::sin(Polygon::vertexTheta[i] + this->getAngle());
		out << "> ,";
	}
	out << "< ";
	out << position[0] + Polygon::vertexR[0] * std::cos(Polygon::vertexTheta[0] + this->getAngle());
	out << ", ";
	out << position[1] + Polygon::vertexR[0] * std::sin(Polygon::vertexTheta[0] + this->getAngle());
	out << ">";
	out << "  texture { finish { ambient 1 diffuse 0 } pigment { color Red} } }" << std::endl;

	return out.str();
}
