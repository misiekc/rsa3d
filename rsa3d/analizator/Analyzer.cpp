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
#include "../../statistics/ArrayFunction.h"
#include "../utils/OMPMacros.h"
#include "../utils/Assertions.h"
#include "DomainAnalyzer.h"

#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>

void Analyzer::Result::print(std::ostream &out) const {
	out << "dir\ttheta\td(theta)\tA\td(A)\td\td(d)\ttheta_inf\td(theta_inf)\tC1\td(C1)\tC2\td(C2)" << std::endl;
	out << dir << "\t" << theta << "\t" << A << "\t" << d << "\t" << thetaInf << "\t" << C1 << "\t" << C2 << std::endl;
}


void Analyzer::fillOrderParameterRange(const RSAShape *s){
	RND rnd(1);
	std::uniform_real_distribution<double> ud(0.0, 2*M_PI);
    auto particle1 = dynamic_cast<const OrderCalculable*>(s);
    std::vector<std::vector<double>> vParameters;
	int total = 1000000;
	_OMP_PARALLEL_FOR
	for(int i=0; i<total; i++){
		RSAShape *other = ShapeFactory::createShape(&rnd);

#if RSA_ANGULAR_DIMENSION==1
		RSAOrientation angle{};
        angle[0] = rnd.nextValue(&ud);
        other->rotate(angle);
#endif
	    auto particle2 = dynamic_cast<const OrderCalculable*>(other);
	    _OMP_CRITICAL(vParameters)
        vParameters.push_back(particle1->calculateOrder(particle2));
	}
	int numberOfParameters = particle1->getNumOfOrderParameters();
	for(int i=0; i<numberOfParameters; i++){
		std::vector<double> statistics;
		statistics.push_back(-std::numeric_limits<double>::infinity());
		statistics.push_back(std::numeric_limits<double>::infinity());
		statistics.push_back(0.0);
		this->orderParameterRange.push_back(statistics);
	}
	for(int i=0; i<total; i++){
		for(int j=0; j<numberOfParameters; j++){
			double d = vParameters[i][j];
			if (this->orderParameterRange[j][0]<d) // max update
				this->orderParameterRange[j][0]=d;
			if (this->orderParameterRange[j][1]>d) // min update
				this->orderParameterRange[j][1]=d;
			this->orderParameterRange[j][2] += d;
		}
	}
	for(int j=0; j<numberOfParameters; j++){
		this->orderParameterRange[j][2] /= total;
	}
}

double Analyzer::applyOrderNormalisation(double d, std::vector<double> &norm){
//	d = 2*(d-norm[2])/(norm[0]-norm[1]); // (d - mean)(max - min)
	return d;
}


void Analyzer::analyzeOrder(const Packing &packing, const NeighbourGrid<const RSAShape> &ng, std::vector<Plot*> *order) {
	if (!this->isOrderCalculable(packing.front()))
    	return;

	double da[RSA_SPATIAL_DIMENSION];
   	double dMax = (*order)[0]->getMax(), dMax2 = dMax*dMax;
   	if (this->orderParameterRange.size()==0)
   		this->fillOrderParameterRange(packing[0]);
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
    		for(std::size_t k{}; k<this->params->packingFractionSurfaceDimension(); k++){
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
    				double d = orderParameters[k];
    				d = this->applyOrderNormalisation(d, this->orderParameterRange[k]);
    				_OMP_CRITICAL(order)
                    (*order)[k]->add(std::sqrt(dist2), d);
    			}
    		}
    	}
    }
}

void Analyzer::analyzeCorrelations(const Packing &packing, const NeighbourGrid<const RSAShape> &ng, Plot *corr){
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
			for (std::size_t k{}; k < this->params->packingFractionSurfaceDimension(); k++){
				da[k] = std::fabs(posi[k] - posj[k]);
				if (da[k]>0.5*this->params->surfaceSize)
					da[k] = this->params->surfaceSize - da[k];
				dist2 += da[k]*da[k];
				if (da[k] > dMax)
					break;
			}
			if (dist2 < dMax2){
				_OMP_CRITICAL(corr)
				corr->add(std::sqrt(dist2));
			}
		}
	}
}

const RSAShape * Analyzer::getNearestNeighbour(const RSAShape *s, const NeighbourGrid<const RSAShape> &ng){
	const RSAShape *nearestNeighbour = nullptr;
	double da[RSA_SPATIAL_DIMENSION];
	double nearestNeighbourDistance2 = std::numeric_limits<double>::max();
	std::vector<const RSAShape *> neighbours;
	RSAVector pos = s->getPosition();
	ng.getNeighbours(&neighbours, pos);
	for(const RSAShape* sn : neighbours){
		if (sn == s)
			continue;
		RSAVector posn = sn->getPosition();
		double dist2 = 0.0;
		for (std::size_t k{}; k<this->params->packingFractionSurfaceDimension(); k++){
			da[k] = std::fabs(pos[k] - posn[k]);
			if (da[k]>0.5*this->params->surfaceSize)
				da[k] = this->params->surfaceSize - da[k];
			dist2 += da[k]*da[k];
		}
		if (nearestNeighbourDistance2 > dist2){
			nearestNeighbourDistance2 = dist2;
			nearestNeighbour = sn;
		}
	}
	return nearestNeighbour;
}

void Analyzer::analyzeNearestNeighbours(const Packing &packing, const NeighbourGrid<const RSAShape> &ng, Plot *minDist){
	double da[RSA_SPATIAL_DIMENSION];
	_OMP_PARALLEL_FOR
   	for(size_t i=0; i<packing.size(); i++){
   		const RSAShape *si = packing[i];
   		const RSAShape *nearestNeighbour = this->getNearestNeighbour(si, ng);
   		if(nearestNeighbour==nullptr)
   			continue;
   		const RSAShape *nearestNeighbourOfNearestNeighbour = this->getNearestNeighbour(nearestNeighbour, ng);

   		if (
   				(nearestNeighbourOfNearestNeighbour==si && si->no < nearestNeighbour->no) ||
				(nearestNeighbourOfNearestNeighbour!=si)
			){
			double dist2 = 0.0;
   			RSAVector posi = si->getPosition();
   			RSAVector posnn = nearestNeighbour->getPosition();
			for(std::size_t k{}; k < this->params->packingFractionSurfaceDimension(); k++){
				da[k] = std::fabs(posi[k] - posnn[k]);
				if (da[k]>0.5*this->params->surfaceSize)
					da[k] = this->params->surfaceSize - da[k];
				dist2 += da[k]*da[k];
			}
			_OMP_CRITICAL(closestNeighbour)
			minDist->add(std::sqrt(dist2));
		}

	}
}
/*
double Analyzer::getPeriodicDistance(const RSAShape *shape1, const RSAShape *shape2) const {
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
/*
void Analyzer::analyzePacking(const Packing &packing, LogPlot *nvt, Plot *asf, double surfaceFactor){
	double lastt = 0.0;
	for(size_t i=0; i<packing.size(); i++){
		double t = packing[i]->time;
		if (nvt != nullptr)
			nvt->addBetween(lastt, t, i);
		if (asf != nullptr){
			double packingFraction = i*surfaceFactor;
			double tries = (t-lastt)/surfaceFactor;
			if (tries<1.0)
				tries = 1.0;
			asf->add(packingFraction, tries);
		}
		lastt = t;
	}
	if (nvt != nullptr)
		nvt->addBetween(lastt, nvt->getMax()+1.0, packing.size());
}
*/
double Analyzer::analyzePacking(const Packing &packing, LogPlot *nvt, Plot *asf, double surfaceFactor){
	double lastt = 0.0;
	double sumOfParticlesVolume = 0.0;
	for(size_t i=0; i<packing.size(); i++){
		double t = packing[i]->time;
		if (nvt != nullptr)
			nvt->addBetween(lastt, t, sumOfParticlesVolume*surfaceFactor);
		if (asf != nullptr){
			double packingFraction = sumOfParticlesVolume*surfaceFactor;
			double tries = (t-lastt)/surfaceFactor;
			if (tries<1.0)
				tries = 1.0;
			asf->add(packingFraction, tries);
		}
		double particleVolume;
		// As for now, we have 2 cases:
		// - surface is flat and can be lower dimensional than particle. Then the first branch correctly counts number
		//   of particles per unit volume assuming that particle volume is 1. In the second branch
		//   params->packingFractionSurfaceDimension() returns just params->surfaceDimension so everything is ok
		// - surface is curved and it cannot be lower dimensional than particle. The first branch counts number of
		//   particles per surface projection volume, and the second the correct projection packing fraction, because
		//   params->packingFractionSurfaceDimension() returns RSA_SPATIAL_DIMENSION - 1
		if (this->params->coverageByNumber)
			particleVolume = packing[i]->getVolume(RSA_SPATIAL_DIMENSION);
		else
			particleVolume = packing[i]->getVolume(this->params->packingFractionSurfaceDimension());
		sumOfParticlesVolume += particleVolume;
		lastt = t;
	}
	if (nvt != nullptr)
		nvt->addBetween(lastt, nvt->getMax()+1.0, sumOfParticlesVolume*surfaceFactor);
	return sumOfParticlesVolume*surfaceFactor;
}


void Analyzer::printKinetics(LogPlot &nvt, std::string filename, double* fixedA, double surfaceFactor, double dTime, Result *res){
	double parameterD = 0.0, errorD = 0.0;
	double **points = new double*[nvt.size()];
	for(int i=0; i<nvt.size(); i++)
		points[i] = new double[2];
	nvt.getAsPoints(points);
	ArrayFunction theta(points, nvt.size());
	int minJ = 0, maxJ = nvt.size()-1;
	double diff, maxX = points[0][maxJ];
	PowerRegression pr;
	LinearRegression lr1;
	LinearRegression lr2;
	std::ofstream file(filename);
 	for (int i = 1, lasti =0; i < nvt.size(); i++) {
		if (points[i][1]==0) // no data
			continue;
		diff = (points[i][1] - points[lasti][1])/(points[i][0] - points[lasti][0]);
		if (diff==0) // no change in data
			continue;
		lasti = i;
	    file.precision(std::numeric_limits<double>::digits10 + 1);
		file << points[i][0] << "\t" << points[i][1] << "\t" << diff << "\t" << theta.get(2.0*points[i][0]) - theta.get(points[i][0]);
		if (points[i][0]>100){ // calculate slope for times in (points[i][0]/100, points[i][0])
			maxJ = i;
			pr.clear();
			lr1.clear();
			lr2.clear();
			maxX = points[maxJ][0];
			int j=minJ, lastj;
			while( !(points[j][0] > 0.01*maxX || (maxX>1.0e+10 && points[j][0] > 0.0000001*maxX)) ){
				j++;
			}
			minJ = j;
			for(j=minJ, lastj=minJ-1; j<maxJ; j++){
				diff = (points[j][1] - points[lastj][1])/(points[j][0] - points[lastj][0]);
				if(diff==0.0)
					continue;
				pr.addXY(points[j][0], diff);
				lastj = j;
			}
			if (pr.size()>=10){
				pr.calculate();
				for(int j=minJ; j<maxJ; j++){
					if (points[j][0]!=0){
						if (fixedA== nullptr){
							lr1.addXY(pow(points[j][0], pr.getA()+1 + pr.getSA()), points[j][1]);
							lr2.addXY(pow(points[j][0], pr.getA()+1 - pr.getSA()), points[j][1]);
						}else{
							lr1.addXY(pow(points[j][0], fixedA[0] + fixedA[1]), points[j][1]);
							lr2.addXY(pow(points[j][0], fixedA[0] - fixedA[1]), points[j][1]);
						}
					}
				}
				lr1.calculate();
				lr2.calculate();
                double A = (lr1.getA()+lr2.getA())/2.0;
                double dA = fabs(A - lr1.getA());
                double B = (lr1.getB()+lr2.getB())/2.0;
                double dB = fabs(B - lr1.getB());
				res->d.value = -1.0 / (pr.getA()+1);
				res->d.error = -pr.getSA()/(pr.getA()+1) * res->d.value;
				if (maxX<dTime){
					parameterD = res->d.value;
					errorD = res->d.error;
				}
				res->thetaInf.value = B;
				res->thetaInf.error = dB;
				res->A.value = -A;
				res->A.error = dA;

				file << "\t" << res->A << "\t" << res->d << "\t" << res->thetaInf << std::endl;
			}else{
				file << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << std::endl;
			}
		}else{
			file << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << std::endl;
		}
	}
	file.close();
	for(int i=0; i<nvt.size(); i++)
		delete[] points[i];
	delete[] points;
	res->d.value = parameterD;
	res->d.error = errorD;
}

void Analyzer::printASF(Plot &asf, std::string filename, int counter, double packingFraction, Result *res){
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
	res->C1.value = asfreg.getC1();
	res->C1.error = asfreg.getSC1();
	res->C2.value = asfreg.getC2();
	res->C2.error = asfreg.getSC2();
}

void Analyzer::printHistogram(Plot& plot, std::string filename){
	double **plotPoints = new double*[plot.size()];
	for(int i=0; i<plot.size(); i++){
		plotPoints[i] = new double[2];
	}
	plot.getAsHistogramPoints(plotPoints);

    std::ofstream file(filename);

    for (int i = 0; i < plot.size()-1; i++) {
		file << plotPoints[i][0] << "\t" << plotPoints[i][1];
		file << std::endl;
	}
	file.close();

	for(int i=0; i<plot.size(); i++){
		delete[] plotPoints[i];
	}
	delete[] plotPoints;
}

void Analyzer::printCorrelations(Plot& correlations, std::string filename){
	double **correlationPoints = new double*[correlations.size()];
	for(int i=0; i<correlations.size(); i++){
		correlationPoints[i] = new double[2];
	}
	correlations.getAsHistogramPoints(correlationPoints);

	double packingVolume;

    std::size_t dim = this->params->packingFractionSurfaceDimension();
	if (dim==2)
		packingVolume = M_PI * pow(correlations.getMax(), 2.0);
	else if (dim==3)
		packingVolume = 4.0/3.0 * M_PI * pow(correlations.getMax(), 3.0);
	else
		return;

    std::ofstream file(filename);

    unsigned long int totalPoints = correlations.getTotalNumberOfHistogramPoints();

    double r{};
	for (int i = 0; i < correlations.size()-1; i++) {
	    double thinShellVolume{};
		if (i>0){
			if (dim==2)
				thinShellVolume = -M_PI * r*r;
			else if (dim==3)
				thinShellVolume = -4.0/3.0* M_PI * r*r*r;
		}
		r = 0.5*(correlationPoints[i][0]+correlationPoints[i+1][0]);
		if (dim==2)
			thinShellVolume += M_PI * r*r;
		else if (dim==3)
			thinShellVolume += 4.0/3.0* M_PI * r*r*r;

		double expectedNumberOfParticles = totalPoints * thinShellVolume / packingVolume;
		double G = correlationPoints[i][1] / expectedNumberOfParticles;
		double dG = std::sqrt(static_cast<double>(correlationPoints[i][1])) / expectedNumberOfParticles;
		
		file << correlationPoints[i][0] << "\t" << G << "\t" << dG;
		file << std::endl;
	}
	file.close();

	for(int i=0; i<correlations.size(); i++){
		delete[] correlationPoints[i];
	}
	delete[] correlationPoints;
}

void Analyzer::printOrder(const std::vector<Plot*> &order, const std::string &filename) const {
    if (order.empty()) return;

	double ***orderPoints = new double**[order.size()];
	for (unsigned short i=0; i<order.size(); i++){
		orderPoints[i] = new double*[order[i]->size()];
		for(unsigned short j=0; j<order[i]->size(); j++){
			orderPoints[i][j] = new double[2];
		}
		order[i]->getAsPoints(orderPoints[i]);
	}

	std::ofstream file(filename);

	for (int i = 0; i < order[0]->size()-1; i++) {
		file << orderPoints[0][i][0];
		for (unsigned short j=0; j<order.size(); j++)
			file << "\t" << orderPoints[j][i][1];
		file << std::endl;
	}
	file.close();

	for(unsigned short i=0; i<order.size(); i++){
		for(int j=0; j<order[i]->size(); j++){
			delete[] orderPoints[i][j];
		}
		delete[] orderPoints[i];
	}
	delete[] orderPoints;
}

void Analyzer::analyzePackingsInDirectory(const std::string &dirName, double mintime, double particleSize,
                                          double correlationsRange) {
    Expects(mintime > 0);
    Expects(particleSize > 0);
    Expects(correlationsRange >= 0);

    auto packingPaths = PackingGenerator::findPackingsInDir(dirName);

    double minmaxTimes[2];
    this->findMinMaxTimes(minmaxTimes, packingPaths);
//    double maxTime = 1.0e+15;

    LogPlot nvt(mintime, minmaxTimes[0], 500);
    Plot asf(0.0, 1.0, 500);
    Plot correlations(0.0, correlationsRange, 500);
    Plot minDist(0.0, correlationsRange, 500);
    std::vector<Plot*> order = this->getFilledOrderVector(correlationsRange);

    double spf = 0.0, spf2 = 0.0;
    int counter = 0;
    double packingSize = this->params->sufraceVolume();

    std::cout << "[Analyzer::analyzePackingsInDirectory] " << std::flush;
	for (const auto &packingPath : packingPaths) {
        Packing packing;
        packing.restore(packingPath);
        double pf = this->analyzePacking(packing, &nvt, &asf, particleSize/packingSize);
        spf += pf;
        spf2 += pf*pf;
        counter++;
        if (correlationsRange>0.0){
            // We want to create a neighbour grid to optimize correlations calculation for "packing fraction" dimension
            // - if some spatial dimensions are projected when calculating the packing fraction, they should also
            // be projected when calculating correlations
            NeighbourGrid<const RSAShape> ng(params->packingFractionSurfaceDimension(), params->surfaceSize, correlationsRange);
            for(const RSAShape *s: packing){
                ng.add(s, s->getPosition());
            }
            this->analyzeCorrelations(packing, ng, &correlations);
            this->analyzeOrder(packing, ng, &order);
            this->analyzeNearestNeighbours(packing, ng, &minDist);
        }
        std::cout << ".";
        std::cout.flush();
	}
	std::cout << std::endl;

	Result result;
	double dTime = (minmaxTimes[0]>0.01*minmaxTimes[1])? minmaxTimes[0] : minmaxTimes[0]/10;
	this->printKinetics(nvt, (dirName + "_kinetics.txt"), nullptr, particleSize / packingSize, dTime, &result);
	double packingFraction = result.thetaInf.value;
	this->printASF(asf, dirName + "_asf.txt", counter, packingFraction, &result);

	result.dir = dirName;
	result.theta.value = spf/counter;
	result.theta.error = sqrt( (spf2/counter - spf*spf/(counter*counter)) / counter);

	result.print(std::cout);

	if (correlationsRange > 0) {
        this->printCorrelations(correlations, dirName + "_cor.txt");
        this->printOrder(order, dirName + "_order.txt");
        this->printHistogram(minDist, dirName + "_dist.txt");
    }
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
        result.push_back(new Plot(0.0, range, 500));
    return result;
}

bool Analyzer::isOrderCalculable(const RSAShape *shape) const {
    return dynamic_cast<const OrderCalculable*>(shape) != nullptr;
}

void Analyzer::findMinMaxTimes(double* minmaxTimes, const std::vector<std::string> &packingPaths) {
    if (this->params->maxTime < std::numeric_limits<double>::infinity()){
        minmaxTimes[0] = this->params->maxTime;
        minmaxTimes[1] = this->params->maxTime;
        return;
    }

    minmaxTimes[1] = 0;
    minmaxTimes[0] = std::numeric_limits<double>::infinity();

    std::cout << "[Analyzer::findMinMaxTime] " << std::flush;
    for (const auto &packingFilename : packingPaths) {
        Packing packing;
        packing.restore(packingFilename);

        double packingTime = packing.back()->time;
        if (packingTime > minmaxTimes[1])
            minmaxTimes[1] = packingTime;
        if (packingTime < minmaxTimes[0])
        	minmaxTimes[0] = packingTime;

        std::cout << "." << std::flush;
    }
    std::cout << " times found: " << minmaxTimes[0] << "\t" << minmaxTimes[1] << std::endl;
}
