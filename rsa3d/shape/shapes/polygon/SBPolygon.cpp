/*
 * SBPolygon.cpp
 *
 *  Created on: 04.06.2018
 *      Author: ciesla
 */

#include "SBPolygon.h"

#include <cmath>

SBPolygon::SBPolygon() : Polygon(){
}

SBPolygon::~SBPolygon() {
}

void SBPolygon::initClass(const std::string &args){
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
	x[2] = x[1] + width*s;	y[2] = y[1] - width*c;
	x[3] = - x[0]; 			y[3] = y[0];
	x[4] = x[2];			y[4] = -y[2];
	x[5] = x[1];			y[5] = -y[1];

	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);
	out << "6 rt ";


	for (size_t i=0; i<6; i++){
		double r, t;
		r = std::sqrt(x[i]*x[i] + y[i]*y[i]);
		t = std::atan2(y[i], x[i]);
		out << r << " " << t << " ";
	}
	out << "5 0 2 1 3 0 3 0 4 3 5";

	Polygon::initClass(out.str());
}
