/*
 * Analyzer.cpp
 *
 *  Created on: 25.06.2017
 *      Author: ciesla
 */

#include "Analyzer.h"
#include "../Shape.h"
#include "../ShapeFactory.h"
#include "../shapes/Cuboid.h"
#include "../../statistics/Plot.h"
#include "../../statistics/LogPlot.h"
#include "../../statistics/LinearRegression.h"
#include "../../statistics/PowerRegression.h"
#include "../../statistics/ASFRegression.h"
#include <fstream>
#include <iostream>
#include <dirent.h>
#include <string.h>
#include <string>

#include <cmath>

Analyzer::Analyzer(Parameters *params) {
	this->params = params;
}

Analyzer::~Analyzer() {
}

std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> * Analyzer::fromFile(std::string filename){
	std::ifstream file(filename, std::ios::binary);
	std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> * v = new std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *>;
	RND rnd(1);

	while(!file.eof()){
		Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *s = ShapeFactory::createShape(&rnd);
		s->restore(file);
		v->push_back(s);
	}
	v->pop_back();

	file.close();
	return v;
}

double P2(double x){
	return 0.5*(3*x*x - 1);
}

double P4(double x){
	return 0.125*(x*x*(35*x*x-30)+3);
}

void Analyzer::calculateOrderParameters(double *result, Cuboid *c1, Cuboid *c2){
	double axisAlignment[3];
	Matrix<3,3> product = (c1->getOrientation().transpose()) * (c2->getOrientation());

	double tmp = 0, row[2], column[2], prod4 = 0;
	axisAlignment[0] = 0.0; row[0] = 0; column[0] = 0;
	for(size_t i = 0; i<3; i++){
		for(size_t j = 0; j<3; j++){
			tmp = fabs(product(i, j));
			prod4 += tmp*tmp*tmp*tmp;
			if (axisAlignment[0] < tmp){
				axisAlignment[0] = tmp;
				row[0] = i;
				column[0] = j;
			}
		}
	}
	axisAlignment[1] = 0.0; row[1] = 0; column[1] = 0;
	for(size_t i = 0; i<3; i++){
		if (i==row[0])
			continue;
		for(size_t j = 0; j<3; j++){
			if (j==column[0])
				continue;
			tmp = fabs(product(i, j));
			if (axisAlignment[1] < tmp){
				axisAlignment[1] = tmp;
				row[1] = i;
				column[1] = j;
			}
		}
	}
	axisAlignment[2] = fabs(product(3 - row[0] - row[1], 3 - column[0] - column[1]));

	result[0] = P2(axisAlignment[0]);
	result[1] = ( P2(axisAlignment[0]) + P2(axisAlignment[1]) + P2(axisAlignment[2]) )/3.0;
	result[2] = P4(axisAlignment[0]);
	result[3] = ( P4(axisAlignment[0]) + P4(axisAlignment[1]) + P4(axisAlignment[2]) )/3.0;
	result[4] = (5.0*prod4 - 9.0)/6.0;
}

void Analyzer::analyzeOrder(std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> *packing, Plot **order){
	double *posi, *posj, da[RSA_SPATIAL_DIMENSION];
	double orderParameters[5];
	for(uint i=0; i<packing->size(); i++){
		posi = (*packing)[i]->getPosition();
		for(uint j=i+1; j<packing->size(); j++){
			posj = (*packing)[j]->getPosition();
			double dist = 0.0;
			for(unsigned short k=0; k<RSA_SPATIAL_DIMENSION; k++){
				da[k] = fabs(posi[k] - posj[k]);
				if (da[k]>0.5*this->params->surfaceSize)
					da[k] = this->params->surfaceSize - da[k];
				dist += da[k]*da[k];
			}
			dist = sqrt(dist);
			if (dist > order[0]->getMax())
				continue;
			this->calculateOrderParameters(orderParameters, (Cuboid *)((*packing)[i]), (Cuboid *)((*packing)[j]) );
			for(ushort k = 0; k<5; k++){
				order[k]->add(dist, orderParameters[k]);
			}
		}
	}
}

void Analyzer::analyzePacking(std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> *packing, LogPlot *nvt, Plot *asf, Plot *corr, double surfaceFactor){
	double lastt = 0.0;
	for(uint i=0; i<packing->size(); i++){
		double t = (*packing)[i]->time;
		if (nvt != NULL)
			nvt->addBetween(lastt, t, i);
		if (asf != NULL){
			double packingFraction = i*surfaceFactor;
			double tries = (t-lastt)/surfaceFactor;
			if (tries<1.0)
				tries = 1.0;
			asf->add(packingFraction, tries);
		}
		lastt = t;
	}
	if (nvt != NULL)
		nvt->addBetween(lastt, nvt->getMax()+1.0, packing->size());
	if (corr!=NULL){
		double *posi, *posj, da[RSA_SPATIAL_DIMENSION];
		for(uint i=0; i<packing->size(); i++){
			posi = (*packing)[i]->getPosition();
			for(uint j=i+1; j<packing->size(); j++){
				posj = (*packing)[j]->getPosition();
				double dist = 0.0;
				for(unsigned short k=0; k<RSA_SPATIAL_DIMENSION; k++){
					da[k] = fabs(posi[k] - posj[k]);
					if (da[k]>0.5*this->params->surfaceSize)
						da[k] = this->params->surfaceSize - da[k];
					dist += da[k]*da[k];
				}
				corr->add(sqrt(dist), 1.0);
			}
		}
	}
}

double * Analyzer::printNvT(LogPlot &nvt, std::string filename, double* fixedA, double surfaceFactor, double *res){
	double **points = new double*[nvt.size()];
	for(int i=0; i<nvt.size(); i++)
		points[i] = new double[2];
	nvt.getAsPoints(points);
	double lastX = points[0][0], lastY = points[0][1], maxX = points[nvt.size()-1][0], diff;
	PowerRegression pr;
	std::ofstream file(filename);
	for (int i = 0; i < nvt.size(); i++) {
		if (points[i][1]!=lastY){
			diff = (points[i][1] - lastY)/(points[i][0] - lastX);
			lastX = points[i][0];
			lastY = points[i][1];
			if (points[i][0] > maxX/100.0){
				pr.addXY(points[i][0], diff);
			}
			file << points[i][0] << "\t" << points[i][1]*surfaceFactor << "\t" << diff*surfaceFactor << std::endl;
		}
	}
	file.close();
	pr.calculate();
	LinearRegression lr1;
	LinearRegression lr2;
	for(int i=0; i<nvt.size(); i++){
		if (points[i][0] > maxX/100.0){
			if (fixedA==NULL){
				lr1.addXY(pow(points[i][0], pr.getA()+1 + pr.getSA()), points[i][1]);
				lr2.addXY(pow(points[i][0], pr.getA()+1 - pr.getSA()), points[i][1]);
			}else{
				lr1.addXY(pow(points[i][0], fixedA[0] + fixedA[1]), points[i][1]);
				lr2.addXY(pow(points[i][0], fixedA[0] - fixedA[1]), points[i][1]);
			}
		}
	}

	for(int i=0; i<nvt.size(); i++)
		delete[] points[i];
	delete[] points;



	lr1.calculate();
	lr2.calculate();
//	ps.println();
	double B = (lr1.getB()+lr2.getB())/2.0;
	double dB = fabs(B - lr1.getB());
	res[0] = pr.getA()+1;
	res[1] = pr.getSA();
	res[2] = B;
	res[3] = dB;
	return res;
}

double * Analyzer::printASF(Plot &asf, std::string filename, int counter, double packingFraction, double *res){
	double **points = new double*[asf.size()];
	for(int i=0; i<asf.size(); i++)
		points[i] = new double[2];
	asf.getAsPoints(points);

	std::ofstream file(filename);
	ASFRegression asfreg;
	for (int i = 0; i < asf.size(); i++) {
		if (points[i][1]>0.0){
			if (points[i][0] < 0.2*packingFraction){
				asfreg.addXY(points[i][0], 1.0 / points[i][1]);
			}
			file << points[i][0] << "\t" << 1.0 / points[i][1] << std::endl;
		}
	}
	file.close();
	for(int i=0; i<asf.size(); i++)
		delete[] points[i];
	delete[] points;
	asfreg.calculate();
	res[0] = asfreg.getC1();
	res[1] = asfreg.getSC1();
	res[2] = asfreg.getC2();
	res[3] = asfreg.getSC2();

	return res;
}

void Analyzer::printCorrelations(Plot& correlations, std::string filename, int counter, double particleSize, double packingFraction){
	double **correlationPoints = new double*[correlations.size()];
	for(int i=0; i<correlations.size(); i++){
		correlationPoints[i] = new double[2];
	}
	correlations.getAsHistogramPoints(correlationPoints);

	double r = 0.0, volume, expectedNumberOfParticles, packingVolume;

	if (RSA_SPATIAL_DIMENSION==2)
		packingVolume = M_PI * pow(correlations.getMax(), 2.0);
	else if (RSA_SPATIAL_DIMENSION==3)
		packingVolume = 4.0/3.0 * M_PI * pow(correlations.getMax(), 3.0);
	else
		return;

	int totalPoints = correlations.getTotalNumberOfPoints();

	std::ofstream file(filename);

	volume = 0.0;
	for (int i = 0; i < correlations.size()-1; i++) {
		if (i>0){
			if (RSA_SPATIAL_DIMENSION==2)
				volume = -M_PI * r*r;
			else if (RSA_SPATIAL_DIMENSION==3)
				volume = -4.0/3.0* M_PI * r*r*r;
		}
		r = 0.5*(correlationPoints[i][0]+correlationPoints[i+1][0]);
		if (RSA_SPATIAL_DIMENSION==2)
			volume += M_PI * r*r;
		else if (RSA_SPATIAL_DIMENSION==3)
			volume += 4.0/3.0* M_PI * r*r*r;

		expectedNumberOfParticles = totalPoints * volume / packingVolume;
		file << correlationPoints[i][0] << "\t" << correlationPoints[i][1] / expectedNumberOfParticles;
		file << std::endl;
	}
	file.close();

	for(int i=0; i<correlations.size(); i++){
		delete[] correlationPoints[i];
	}
	delete[] correlationPoints;
}
void Analyzer::printOrder(Plot **order, std::string filename){
	double **orderPoints[5];
	for (ushort i=0; i<5; i++){
		orderPoints[i] = new double*[order[i]->size()];
		for(ushort j=0; j<order[i]->size(); j++){
			orderPoints[i][j] = new double[2];
		}
		order[i]->getAsPoints(orderPoints[i]);
	}

	std::ofstream file(filename);

	for (int i = 0; i < order[0]->size()-1; i++) {
		file << orderPoints[0][i][0];
		for (ushort j=0; j<5; j++)
			file << "\t" << orderPoints[j][i][1];
		file << std::endl;
	}
	file.close();

	for(ushort i=0; i<5; i++){
		for(int j=0; j<order[i]->size(); j++){
			delete[] orderPoints[i][j];
		}
		delete[] orderPoints[i];
	}
}

void Analyzer::analyzePackingsInDirectory(char *sdir, double mintime, double particleSize){
	dirent *de;
	char prefix[] = "packing";
	LogPlot *nvt = new LogPlot(mintime, this->params->maxTime, 200);
	Plot *asf = new Plot(0.0, 0.6, 200);
	Plot *correlations = new Plot(0.0, 10.0, 200);
	Plot *order[5];
	for(ushort i=0; i<5; i++)
		order[i] = new Plot(0.0, 10.0, 200);

	double n = 0.0, n2 = 0.0;
	int counter = 0;
	double packingSize = pow(params->surfaceSize, RSA_SPATIAL_DIMENSION);
	DIR *dir = opendir(sdir);
	std::string dirname(sdir);
	while ((de = readdir(dir)) != NULL){
		if (strncmp(prefix, de->d_name, strlen(prefix))==0){
			std::string filename(de->d_name);
 			std::vector<Shape<RSA_SPATIAL_DIMENSION, RSA_ANGULAR_DIMENSION> *> * packing = this->fromFile(dirname + "/" + filename);
			n += packing->size();
			n2 += packing->size()*packing->size();
			counter++;
			this->analyzePacking(packing, nvt, asf, correlations, particleSize/packingSize);
			if (this->params->particleType.compare("Cuboid")==0)
				this->analyzeOrder(packing, order);
			delete packing;
			std::cout << ".";
			std::cout.flush();
		}
	}
	(void)closedir(dir);
	std::cout << std::endl;
	double res[4];

	//	double fixedA[] = {1.0, 0.0};

	this->printNvT(*nvt, (dirname + "_nvt.txt"), NULL, particleSize/packingSize, res);
	std::cout << sdir << "\t" << n/counter << "\t" << sqrt( (n2/counter - n*n/(counter*counter)) / counter) << "\t" <<
			 	 res[0] << "\t" << res[1] << "\t" << res[2] << "\t" << res[3] << "\t";
	delete nvt;
	double packingFraction = res[2]*particleSize / packingSize;
	this->printASF(*asf, dirname + "_asf.txt", counter, packingFraction, res);

	std::cout << res[0] << "\t" << res[1] << "\t" << res[2] << "\t" << res[3] <<  std::endl;

	delete asf;

	this->printCorrelations(*correlations, dirname + "_cor.txt", counter, particleSize, packingFraction);
	delete correlations;

	this->printOrder(order, dirname + "_order.txt");
	for(ushort i=0; i<5; i++)
		delete order[i];
}
