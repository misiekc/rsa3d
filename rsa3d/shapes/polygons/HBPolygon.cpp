/*
 * HBPolygon.cpp
 *
 *  Created on: 04.06.2018
 *      Author: ciesla
 */

#include "HBPolygon.h"

HBPolygon::HBPolygon() : Polygon(){
}

HBPolygon::~HBPolygon() {
}

void HBPolygon::initClass(const std::string &args){
	std::istringstream in(args);

	double width, alpha;
	in >> width;
	in >> alpha;

	double s = sin(alpha);
	double c = cos(alpha);
	double d = width / s;

	double x[6], y[6];

	x[0] = -0.5*d;			y[0] = 0.0;
	x[1] = x[0] + c;		y[1] = y[0] + s;
	x[2] = x[1] + d;		y[2] = y[1];
	x[3] = - x[0]; 			y[3] = y[0];
	x[4] = x[2];			y[4] = -y[2];
	x[5] = x[1];			y[5] = -y[1];

	for (size_t i=0; i<6; i++){
		double r, t;
		r = std::sqrt(x[i]*x[i] + y[i]*y[i]);
		t = std::atan2(y[i], x[i]);
		Polygon::vertexR.push_back(r);
		Polygon::vertexTheta.push_back(t);
	}
	std::pair<size_t, size_t> s0(0, 2);
	std::pair<size_t, size_t> s1(1, 3);
	std::pair<size_t, size_t> s2(0, 3);
	std::pair<size_t, size_t> s3(0, 4);
	std::pair<size_t, size_t> s4(3, 5);

	Polygon::helperSegments.push_back(s0);
	Polygon::helperSegments.push_back(s1);
	Polygon::helperSegments.push_back(s2);
	Polygon::helperSegments.push_back(s3);
	Polygon::helperSegments.push_back(s4);

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
