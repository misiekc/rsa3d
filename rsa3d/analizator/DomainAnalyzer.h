//
// Created by Michal Ciesla on 11/14/22.
//

#ifndef RSA3D_DOMAINANALYZER_H
#define RSA3D_DOMAINANALYZER_H

#include <unordered_set>
#include "../Packing.h"
#include "../shape/shapes/Rectangle.h"
#include "../shape/shapes/DiscreteOrientationsShape2_1.h"
#include "../../statistics/Plot.h"
#include "../Parameters.h"
#include "../NeighbourGrid.h"
#include "../boundary_conditions/PeriodicBC.h"
#include "../../statistics/LogPlot.h"

class Domain {
private:
    NeighbourGrid<const RSAShape> *ng;
    RSABoundaryConditions *bc;
    double orientation{};
    double boxLength;
    double rectangleSides[2] = {Rectangle::longer, Rectangle::shorter};

    // iAxis = 0 - returns shape at the left side of the packing
    // iAxis = 1 - returns shape at the bottom side of the packing
    const RSAShape *getBeginningShape(size_t iAxis) {
        for (const RSAShape *s: this->shapes) {
            if ( (this->orientation == DiscreteOrientationsShape2_1::allowedOrientations[0][0]) &&
                    (std::abs(s->getPosition()[iAxis]) < 0.5 * (rectangleSides[iAxis] + rectangleSides[1])))
                return s;
            else if(
                    (this->orientation != DiscreteOrientationsShape2_1::allowedOrientations[iAxis][0]) &&
                    (std::abs(s->getPosition()[iAxis]) < 0.5 * (rectangleSides[1] + rectangleSides[1])))
                return s;
        }
        return nullptr;
    }

    bool percolationLeftRight(){
        const RSAShape *seed = this->getBeginningShape(0);
        std::vector<const RSAShape *> neighbours;
//        std::unordered_set
//        this->ng->getNeighbours(&neighbours, seed->getPosition());
//        for(const RSAShape *s: neighbours){

//        }
    }


public:
    std::unordered_set<const RSAShape *> shapes;

    explicit Domain(Parameters &params){
        this->bc = new RSAPeriodicBC(params.surfaceSize);
        this->ng = new NeighbourGrid<const RSAShape>(params.surfaceDimension, params.surfaceSize, RSAShape::getNeighbourListCellSize());
        this->boxLength = params.surfaceSize;
    }

    ~Domain(){
        delete this->bc;
        delete this->ng;
    }

    bool testShape(const RSAShape *shape) const{
        if (this->shapes.empty())
            return true;
        if (this->getOrientation() != shape->getOrientation()[0])
            return false;

        std::vector<const RSAShape *> neighbours;
        this->ng->getNeighbours(&neighbours, shape->getPosition());
        for (const RSAShape *s : neighbours){ // we check only neighbours of shapes inside a domain
            RSAVector trans = this->bc->getTranslation(shape->getPosition(), s->getPosition());
            RSAVector dpos = shape->getPosition() - s->getPosition() - trans;
            if (shape->getOrientation()[0] == 0) {
                if (
                        (std::abs(dpos[0]) < Rectangle::longer + Rectangle::shorter)
                        &&
                        (std::abs(dpos[1]) < 2*Rectangle::shorter)
                        )
                    return true;
            } else {
                if (
                        (std::abs(dpos[1]) < Rectangle::longer + Rectangle::shorter)
                        &&
                        (std::abs(dpos[0]) < 2*Rectangle::shorter)
                        )
                    return true;
            }
        }
        return false;
    }

    bool isPercolating() const{
        if (static_cast<double>(this->shapes.size())*(Rectangle::longer + Rectangle::shorter) < this->boxLength)
            return false;
//        return (this->percolateHorizontally() && this->percolateVertically());
    }

    size_t size() const{
        return this->shapes.size();
    }

    double getOrientation() const{
        return this->orientation;
    }

    void addShape(const RSAShape *shape) {
        if (this->shapes.empty())
            this->orientation = shape->getOrientation()[0];
        this->shapes.insert(shape);
        this->ng->add(shape, shape->getPosition());
    }
};

class DomainAnalyzer {
private:
    static Parameters params;
    static double analyzeOrder(const Packing &packing);
    static void addNeighboursToQueue(std::vector<const RSAShape*> &neighbours, std::unordered_set<const RSAShape *> &queue);
    static void dividePackingIntoDomains(const Packing &packing, std::vector<const Domain *> &domains);
    static void printHistogram(Plot &plot, const std::string &filename);
    static void toWolfram(const std::vector<const Domain *> &domains, const std::string &filename);

public:
    static void init(const Parameters &params);
    static void analyzeOrderDirectory(const std::string &dirName);
    static void analyzeDomains(const std::string &dirName);
    static void toWolfram(const Packing& packing, const std::string &filename);
};

#endif //RSA3D_DOMAINANALYZER_H
