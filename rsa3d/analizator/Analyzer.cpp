/*
 * Analyzer.cpp
 *
 *  Created on: 25.06.2017
 *      Author: ciesla
 */

#include "Analyzer.h"
#include "../Shape.h"
#include "../ShapeFactory.h"
#include <Plot.h>
#include <LogPlot.h>
#include <LinearRegression.h>
#include <PowerRegression.h>
#include <ASFRegression.h>
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

std::vector<Shape *> * Analyzer::fromFile(std::string filename){
	std::ifstream file(filename, std::ios::binary);
	std::vector<Shape *> * v = new std::vector<Shape *>;
	RND rnd(1);

	while(!file.eof()){
		Shape *s = ShapeFactory::createShape(&rnd);
		s->restore(file);
		v->push_back(s);
	}
	v->pop_back();

	file.close();
	return v;
}

void Analyzer::analyzePacking(std::vector<Shape *> *packing, LogPlot *nvt, Plot *asf, Plot *corr, double surfaceFactor){
	double lastt = 0.0;
	for(unsigned int i=0; i<packing->size(); i++){
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
		double *posi, *posj, *da = new double[this->params->dimension];
		for(unsigned int i=0; i<packing->size(); i++){
			posi = (*packing)[i]->getPosition();
			for(unsigned int j=i+1; j<packing->size(); j++){
				posj = (*packing)[j]->getPosition();
				double dist = 0.0;
				for(unsigned char k=0; k<this->params->dimension; k++){
					da[k] = fabs(posi[k] - posj[k]);
					if (da[k]>0.5*this->params->surfaceSize)
						da[k] = this->params->surfaceSize - da[k];
					dist += da[k]*da[k];
				}
				corr->add(sqrt(dist), 1.0);
			}
		}
		delete[] da;
	}
}

double * Analyzer::printNvT(LogPlot &nvt, std::string filename, double surfaceFactor, double *res){
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
		if (points[i][0] > maxX/1000.0){
			lr1.addXY(pow(points[i][0], pr.getA()+1+pr.getSA()), points[i][1]);
			lr2.addXY(pow(points[i][0], pr.getA()+1-pr.getSA()), points[i][1]);
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

void Analyzer::printCorr(Plot& corr, std::string filename, int counter, double particleSize, double packingFraction){
	double **points = new double*[corr.size()];
	for(int i=0; i<corr.size(); i++)
		points[i] = new double[2];
	corr.getAsHistogramPoints(points);

	double r = 0.0, volume, expectedNumberOfParticles, packingVolume;

	if (this->params->dimension==2)
		packingVolume = M_PI * pow(corr.getMax(), 2.0);
	else if (this->params->dimension==3)
		packingVolume = 4.0/3.0 * M_PI * pow(corr.getMax(), 3.0);
	else
		return;

	int totalPoints = corr.getTotalNumberOfPoints();

	std::ofstream file(filename);

	volume = 0.0;
	for (int i = 0; i < corr.size()-1; i++) {
		if (i>0){
			if (this->params->dimension==2)
				volume = -M_PI * r*r;
			else if (this->params->dimension==3)
				volume = -4.0/3.0* M_PI * r*r*r;
		}
		r = 0.5*(points[i][0]+points[i+1][0]);
		if (this->params->dimension==2)
			volume += M_PI * r*r;
		else if (this->params->dimension==3)
			volume += 4.0/3.0* M_PI * r*r*r;

		expectedNumberOfParticles = totalPoints * volume / packingVolume;
		file << points[i][0] << "\t" << points[i][1] / expectedNumberOfParticles << std::endl;
	}
	file.close();
	for(int i=0; i<corr.size(); i++)
		delete[] points[i];
	delete[] points;
}

void Analyzer::analyzePackingsInDirectory(char *sdir, double mintime, double particleSize){
	dirent *de;
	char prefix[] = "packing";
	LogPlot *nvt = new LogPlot(mintime, this->params->maxTime, 200);
	Plot *asf = new Plot(0.0, 0.6, 200);
	Plot *corr = new Plot(0.0, 10.0, 200);
	double n = 0.0, n2 = 0.0;
	int counter = 0;
	double packingSize = pow(params->surfaceSize, params->dimension);
	DIR *dir = opendir(sdir);
	std::string dirname(sdir);
	while ((de = readdir(dir)) != NULL){
		if (strncmp(prefix, de->d_name, strlen(prefix))==0){
			std::string filename(de->d_name);
 			std::vector<Shape *> * packing = this->fromFile(dirname + "/" + filename);
			n += packing->size();
			n2 += packing->size()*packing->size();
			counter++;
			this->analyzePacking(packing, nvt, asf, corr, particleSize/packingSize);
			delete packing;
			std::cout << ".";
			std::cout.flush();
		}
	}
	(void)closedir(dir);
	std::cout << std::endl;
	double res[4];

	this->printNvT(*nvt, (dirname + "_nvt.txt"), particleSize/packingSize, res);
	std::cout << sdir << "\t" << n/counter << "\t" << sqrt( (n2/counter - n*n/(counter*counter)) / counter) << "\t" <<
			 	 res[0] << "\t" << res[1] << "\t" << res[2] << "\t" << res[3] << "\t";
	delete nvt;
	double packingFraction = res[2]*particleSize / packingSize;
	this->printASF(*asf, dirname + "_asf.txt", counter, packingFraction, res);

	std::cout << res[0] << "\t" << res[1] << "\t" << res[2] << "\t" << res[3] <<  std::endl;

	delete asf;

	this->printCorr(*corr, dirname + "_cor.txt", counter, particleSize, packingFraction);
	delete corr;
}
