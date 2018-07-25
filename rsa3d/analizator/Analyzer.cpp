/*
 * Analyzer.cpp
 *
 *  Created on: 25.06.2017
 *      Author: ciesla
 */

#include "Analyzer.h"
#include "../shape/ShapeFactory.h"
#include "../../statistics/LinearRegression.h"
#include "../../statistics/PowerRegression.h"
#include "../../statistics/ASFRegression.h"
#include <fstream>
#include <iostream>
#include <dirent.h>
#include <cstring>
#include <cmath>


void Analyzer::analyzeOrder(const Packing &packing, NeighbourGrid<const RSAShape> &ng, std::vector<Plot*> *order) {
    if (!this->isOrderCalculable(packing.front()))
    	return;

  	double da[RSA_SPATIAL_DIMENSION];
   	double dMax = (*order)[0]->getMax(), dMax2 = dMax*dMax;
   	_OMP_PARALLEL_FOR
   	for(size_t i=0; i<packing.size(); i++){
   		const RSAShape *si = packing[i];
    	RSAVector posi = si->getPosition();
   	   	std::vector<const RSAShape *> neighbours;
    	ng.getNeighbours(&neighbours, posi);
    	for(const RSAShape *sj : neighbours){
            if (sj->no <= si->no)
                continue;

    		RSAVector posj = sj->getPosition();
    		double dist2 = 0.0;
    		for(unsigned short k=0; k<RSA_SPATIAL_DIMENSION; k++){
    			da[k] = fabs(posi[k] - posj[k]);
    			if (da[k]>0.5*this->params->surfaceSize)
    				da[k] = this->params->surfaceSize - da[k];
    			dist2 += da[k]*da[k];
    			if (da[k] > dMax)
    				break;
    		}
    		if (dist2 < dMax2){
                auto particle1 = dynamic_cast<const OrderCalculable*>(si);
                auto particle2 = dynamic_cast<const OrderCalculable*>(sj);
                std::vector<double> orderParameters = particle1->calculateOrder(particle2);
    			for(std::size_t k = 0; k < orderParameters.size(); k++){
    				_OMP_CRITICAL(order)
                    (*order)[k]->add(std::sqrt(dist2), orderParameters[k]);
    			}
    		}
    	}
    }
}

void Analyzer::analyzeCorrelations(const Packing &packing, NeighbourGrid<const RSAShape> &ng, Plot *corr){
	double da[RSA_SPATIAL_DIMENSION];
	double dMax = corr->getMax(), dMax2 = dMax*dMax;
	_OMP_PARALLEL_FOR
   	for(size_t i=0; i<packing.size(); i++){
   		const RSAShape *si = packing[i];
    	RSAVector posi = si->getPosition();
       	std::vector<const RSAShape *> neighbours;
    	ng.getNeighbours(&neighbours, posi);
		for(const RSAShape* sj : neighbours){
			if (sj->no <= si->no)
                continue;

			RSAVector posj = sj->getPosition();
			double dist2 = 0.0;
			for(unsigned short k=0; k<RSA_SPATIAL_DIMENSION; k++){
				da[k] = fabs(posi[k] - posj[k]);
				if (da[k]>0.5*this->params->surfaceSize)
					da[k] = this->params->surfaceSize - da[k];
				dist2 += da[k]*da[k];
				if (da[k] > dMax)
					break;
			}
			if (dist2 < dMax2){
				_OMP_CRITICAL(correlations)
				corr->add(std::sqrt(dist2));
			}
		}
	}
}
/*
double Analyzer::getPetiodicDistance(const RSAShape *shape1, const RSAShape *shape2) const {
    RSAVector delta = shape1->getPosition() - shape2->getPosition();
    RSAVector periodicDelta;
    for(unsigned short k=0; k<RSA_SPATIAL_DIMENSION; k++) {
        periodicDelta[k] = fabs(delta[k]);
        if (periodicDelta[k]>0.5*this->params->surfaceSize)
            periodicDelta[k] = this->params->surfaceSize - periodicDelta[k];
    }
    return periodicDelta.norm();
}
*/

void Analyzer::analyzePacking(const Packing &packing, LogPlot *nvt, Plot *asf, double surfaceFactor){
	double lastt = 0.0;
	for(size_t i=0; i<packing.size(); i++){
		double t = packing[i]->time;
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
		nvt->addBetween(lastt, nvt->getMax()+1.0, packing.size());
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

void Analyzer::printCorrelations(Plot& correlations, std::string filename){
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

	unsigned long int totalPoints = correlations.getTotalNumberOfHistogramPoints();

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

void Analyzer::printOrder(const std::vector<Plot*> &order, const std::string &filename) const {
	double ***orderPoints = new double**[order.size()];
	for (ushort i=0; i<order.size(); i++){
		orderPoints[i] = new double*[order[i]->size()];
		for(ushort j=0; j<order[i]->size(); j++){
			orderPoints[i][j] = new double[2];
		}
		order[i]->getAsPoints(orderPoints[i]);
	}

	std::ofstream file(filename);

	for (int i = 0; i < order[0]->size()-1; i++) {
		file << orderPoints[0][i][0];
		for (ushort j=0; j<order.size(); j++)
			file << "\t" << orderPoints[j][i][1];
		file << std::endl;
	}
	file.close();

	for(ushort i=0; i<order.size(); i++){
		for(int j=0; j<order[i]->size(); j++){
			delete[] orderPoints[i][j];
		}
		delete[] orderPoints[i];
	}
	delete[] orderPoints;
}

void Analyzer::analyzePackingsInDirectory(char *sdir, double mintime, double particleSize, double correlationsRange) {
	dirent *de;
	char prefix[] = "packing";
	LogPlot *nvt = new LogPlot(mintime, this->params->maxTime, 200);
	Plot *asf = new Plot(0.0, 0.6, 200);
	Plot *correlations = new Plot(0.0, correlationsRange, 200);
	std::vector<Plot*> order = this->getFilledOrderVector(correlationsRange);

	double n = 0.0, n2 = 0.0;
	int counter = 0;
	double packingSize = pow(params->surfaceSize, RSA_SPATIAL_DIMENSION);
	DIR *dir = opendir(sdir);
	std::string dirname(sdir);
	while ((de = readdir(dir)) != NULL){
		if (strncmp(prefix, de->d_name, strlen(prefix))==0){
			std::string filename(de->d_name);
			Packing packing;
			packing.restore(dirname + "/" + filename);
			n += packing.size();
			n2 += packing.size()*packing.size();
			counter++;
			this->analyzePacking(packing, nvt, asf, particleSize/packingSize);
			NeighbourGrid<const RSAShape> ng(params->surfaceSize, correlationsRange);
			for(const RSAShape *s: packing){
				ng.add(s, s->getPosition());
			}
			this->analyzeCorrelations(packing, ng, correlations);
			this->analyzeOrder(packing, ng, &order);
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

	this->printCorrelations(*correlations, dirname + "_cor.txt");
	delete correlations;

	this->printOrder(order, dirname + "_order.txt");
	for (Plot *orderPlot : order)
		delete orderPlot;
}

std::vector<Plot*> Analyzer::getFilledOrderVector(double range) const {
    RND rnd;
    std::unique_ptr<RSAShape> shape(ShapeFactory::createShape(&rnd));
    if (!this->isOrderCalculable(shape.get()))
        return {};

    std::vector<Plot*> result;
    auto orderCalculable = dynamic_cast<OrderCalculable*>(shape.get());
    std::size_t paramNum = orderCalculable->getNumOfOrderParameters();
    for (std::size_t i = 0; i < paramNum; i++)
        result.push_back(new Plot(0.0, range, 200));
    return result;
}

bool Analyzer::isOrderCalculable(const RSAShape *shape) const {
    return dynamic_cast<const OrderCalculable*>(shape) != nullptr;
}
