//
// Created by pkua on 05.09.2019.
//

#include <fstream>

#include "ContactFunctionMode.h"
#include "../shape/ShapeFactory.h"

void ContactFunctionMode::initializeForArguments(const ProgramArguments &arguments) {
    this->params = arguments.getParameters();
    std::vector<std::string> positionalArguments = arguments.getPositionalArguments();
    if (positionalArguments.size() != 1)
        this->stepAngle=50;
    else
    	this->stepAngle = std::stoi(positionalArguments[0]);
}

void ContactFunctionMode::absoluteTranslate(RSAShape *s, RSAVector *v){
	RSAVector vTmp = s->getPosition();
	s->translate(-vTmp);
	s->translate(*v);
}

double ContactFunctionMode::calculate(RSAShape *s1, RSAShape *s2) {

	const double factor = 1.5;
	const double precision = 1.0e-10;

	RSAVector vPos = s1->getPosition(), vTmp;
	s1->translate(-vPos);
	s2->translate(-vPos);
	vPos = s2->getPosition();

	double dmin = 1.0, dmax = 1.0;
	while(s1->overlap(&(this->bc), s2)){
		dmax *=factor;
		RSAVector vTmp = vPos*dmax;
		this->absoluteTranslate(s2, &vTmp);
	}

	this->absoluteTranslate(s2, &vPos);
	// previous position of s2 restored
	while(!s1->overlap(&(this->bc), s2)){
		dmin /=factor;
		RSAVector vTmp = vPos*dmin;
		this->absoluteTranslate(s2, &vTmp);
	}

	while(dmax-dmin>precision){
		double d = 0.5*(dmin+dmax);
		RSAVector vTmp = vPos*d;
		this->absoluteTranslate(s2, &vTmp);
		if (s1->overlap(&(this->bc), s2)){
			dmin = d;
		}else{
			dmax = d;
		}
	}
	return 0.5*(dmin+dmax);
}

void ContactFunctionMode::calculate2D(std::string filename){

	RND rnd(0);
	RSAShape *s1 = ShapeFactory::createShape(&rnd);
	RSAShape *s2 = ShapeFactory::createShape(&rnd);
	RSAVector s2Pos;
	s2Pos[0] = 1.0; s2Pos[1] = 0.0;
	s2->translate(s2Pos);

	RSAOrientation dAngle{s1->getAngularVoxelSize()/this->stepAngle};

	std::ofstream file(filename);
	for (double a1 = 0; a1<s1->getAngularVoxelSize(); a1+=dAngle[0]){
		for (double a2 = 0; a2<s2->getAngularVoxelSize(); a2+=dAngle[0]){
			double d = calculate(s1, s2);
			file << a1 << "\t" << a2 << "\t" << d << std::endl;
			s2->rotate(dAngle);
		}
		s1->rotate(dAngle);
	}
	file.close();
}


void ContactFunctionMode::run() {
	this->calculate2D(this->params.getPackingSignature() + "_cf.txt");
}

void ContactFunctionMode::printHelp(std::ostream &out, const ProgramArguments &arguments) {
    out << arguments.formatUsage("[precision]") << std::endl;
    out << std::endl;
    out << "It generates contact function data for two-dimensional shapes. See for" << std::endl;
    out << "example '" << arguments.getCmd() << " help simulate' for information about the *.dat file. Note, that it" << std::endl;
    out << "requires to be invoked with the same parameters as were used to generate the packings." << std::endl;
}
