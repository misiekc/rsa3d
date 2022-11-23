//
// Created by Michal Ciesla on 11/14/22.
//

#ifndef RSA3D_DOMAINANALYZER_H
#define RSA3D_DOMAINANALYZER_H

#include <unordered_set>
#include "../Packing.h"
#include "../shape/shapes/Rectangle.h"
#include "../../statistics/Plot.h"
#include "../Parameters.h"
#include "../NeighbourGrid.h"
#include "../boundary_conditions/PeriodicBC.h"

class Domain {
private:
    NeighbourGrid<const RSAShape> *ng;
    BoundaryConditions<RSA_SPATIAL_DIMENSION> *bc;
    double orientation;

public:
    std::unordered_set<const RSAShape *> shapes;

    Domain(Parameters &params){
        this->bc = new PeriodicBC<RSA_SPATIAL_DIMENSION>(params.surfaceSize);
        this->ng = new NeighbourGrid<const RSAShape>(params.surfaceDimension, params.surfaceSize, RSAShape::getNeighbourListCellSize());
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

    size_t size() const{
        return this->shapes.size();
    }

    double getOrientation() const{
        return this->orientation;
    }

    void addShape(const RSAShape *shape) {
        if (this->shapes.size()==0)
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

public:
    static void init(const Parameters &params);
    static void analyzeOrderDirectory(const std::string &dirName);
    static void analyzeDomains(const std::string &dirName);
    static void toWolfram(Packing packing, const std::string &filename);
};

#endif //RSA3D_DOMAINANALYZER_H
