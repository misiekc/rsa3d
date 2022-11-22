//
// Created by ciesla on 11/14/22.
//

#ifndef RSA3D_DOMAINANALYZER_H
#define RSA3D_DOMAINANALYZER_H

#include "../Packing.h"
#include "../shape/shapes/Rectangle.h"
#include "../../statistics/Plot.h"

class Domain {
public:
    std::vector<const RSAShape *> shapes;

    bool testShape(const RSAShape *shape) {
        if (this->shapes.empty())
            return true;
        if (this->shapes[0]->getOrientation()[0] != shape->getOrientation()[0])
            return false;

        for (const RSAShape *s: this->shapes) {
            RSAVector dpos = shape->getPosition() - s->getPosition();
            if (shape->getOrientation()[0] == 0) {
                if (
                        (std::abs(dpos[0]) - Rectangle::longer < Rectangle::shorter)
                        &&
                        (std::abs(dpos[1]) - Rectangle::shorter < Rectangle::shorter)
                        )
                    return true;
            } else {
                if (
                        (std::abs(dpos[1]) - Rectangle::longer < Rectangle::shorter)
                        &&
                        (std::abs(dpos[0]) - Rectangle::shorter < Rectangle::shorter)
                        )
                    return true;
            }
        }
        return false;
    }

    int size() {
        return this->shapes.size();
    }

    double orientation() {
        return this->shapes[0]->getOrientation()[0];
    }

    void addShape(const RSAShape *shape) {
        this->shapes.push_back(shape);
    }

    std::string toWolfram(){
        std::string sRes = "";
        for (const RSAShape *s: this->shapes){
            sRes.append(s->toWolfram());
            sRes.append(",\n");
        }
        return sRes.substr(0, sRes.size()-2);
    }
};

class DomainAnalyzer {
private:
    static double analyzeOrder(const Packing &packing);
    static void dividePackingIntoDomains(const Packing &packing, std::vector<Domain> &domains);
    static void printHistogram(Plot &plot, const std::string &filename);

public:
    static void analyzeOrderDirectory(const std::string &dirName);
    static void exportPacking(const std::string &fileIn, const std::string &fileOut);
    static void analyzeDomains(const std::string &dirName);

};


#endif //RSA3D_DOMAINANALYZER_H
