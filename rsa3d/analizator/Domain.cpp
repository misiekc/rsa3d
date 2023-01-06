//
// Created by ciesla on 12/27/22.
//

#include "Domain.h"

void Domain::addShape(const RSAShape *shape) {
    if(!this->active){
        throw std::runtime_error("cannot modify passive domain");
    }
    this->shapes.push_back(shape);
    this->ng->add(shape, shape->getPosition());
}

Domain::Domain(double surfaceSize, unsigned short int surfaceDimension, const RSAShape *seed){
    this->surfaceSize = surfaceSize;
    this->surfaceDimension = surfaceDimension;
    this->ng = new NeighbourGrid<const RSAShape>(this->surfaceDimension, this->surfaceSize, RSAShape::getNeighbourListCellSize());
    this->orientation = seed->getOrientation()[0];
    this->addShape(seed);
}

Domain::~Domain(){
    if (this->ng != nullptr)
        delete this->ng;
}

bool Domain::testAndAddShape(const RSAShape *shape, RSABoundaryConditions *bc){
    if(!this->active){
        throw std::runtime_error("cannot modify passive domain");
    }
    if (this->orientation != shape->getOrientation()[0])
        return false;

    std::vector<const RSAShape *> neighbours;
    this->ng->getNeighbours(&neighbours, shape->getPosition());
    for (const RSAShape *s : neighbours){ // we check only neighbours of shapes inside a domain
        RSAVector trans = bc->getTranslation(shape->getPosition(), s->getPosition());
        RSAVector dpos = shape->getPosition() - s->getPosition() - trans;
        if (shape->getOrientation()[0] == 0) {
            if (
                    (std::abs(dpos[0]) < Rectangle::longer + Rectangle::shorter)
                    && (std::abs(dpos[1]) < 2 * Rectangle::shorter)
                    ) {
                this->addShape(shape);
                return true;
            }
        } else {
            if (
                    (std::abs(dpos[1]) < Rectangle::longer + Rectangle::shorter)
                    &&
                    (std::abs(dpos[0]) < 2 * Rectangle::shorter)
                    ) {
                this->addShape(shape);
                return true;
            }
        }
    }
    return false;
}

size_t Domain::size() const{
    return this->shapes.size();
}

double Domain::getOrientation() const{
    return this->orientation;
}

std::vector<const RSAShape *> Domain::getShapes() const{
    return this->shapes;
}

void Domain::makePasive() {
    this->active = false;
    delete this->ng;
    this->ng = nullptr;
}

