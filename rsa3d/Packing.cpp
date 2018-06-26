//
// Created by PKua on 17.06.18.
//

#include "Packing.h"
#include "shape/ShapeFactory.h"
#include "Utils.h"
#include <fstream>

Packing::~Packing() {
    for (auto shape : this->packing)
        delete shape;
}

Packing::Packing(const Packing &other) {
    for (auto &shape : other.packing)
        this->addShape(shape->clone());
}

Packing &Packing::operator=(const Packing &other) {
    this->clear();
    for (auto &shape : other.packing)
        this->addShape(shape->clone());
    return *this;
}

Packing::Packing(Packing &&other) noexcept : packing(std::move(other.packing)) {
    other.packing.clear();
}

Packing &Packing::operator=(Packing &&other) noexcept {
    this->clear();
    this->packing = std::move(other.packing);
    other.packing.clear();
    return *this;
}

void Packing::store(std::ostream &out) const {
    for(auto &shape : this->packing)
        shape->store(out);
}

void Packing::restore(std::istream &in) {
    RND rnd(1);
    packing.clear();

    while (!in.eof()) {
        RSAShape *shape = ShapeFactory::createShape(&rnd);
        shape->restore(in);
        this->addShape(shape);
    }
    this->removeShape(this->size() - 1);
}

void Packing::store(const std::string &filename) const {
    std::ofstream file(filename, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot open file " + filename + " to store packing");
    this->store(file);
    file.close();
}

void Packing::restore(const std::string &filename) {
    std::ifstream file(filename, std::ios_base::binary);
    if (!file) throw std::runtime_error("Cannot open file " + filename + " to restore packing");
    this->restore(file);
    file.close();
}

void Packing::expandOnPBC(double linearSize, double expandMargin) {
    for (std::size_t i = 0; i < RSA_SPATIAL_DIMENSION; i++) {
        std::size_t oldSize = this->size();
        for (std::size_t j = 0; j < oldSize; j++) {
            auto shape = (*this)[j];
            RSAVector position = shape->getPosition();
            if (position[i] < expandMargin * linearSize)
                expandShapeOnBC(shape, linearSize, i);
            else if (position[i] > (1 - expandMargin) * linearSize)
                expandShapeOnBC(shape, -linearSize, i);
        }
    }
}

/* Helper method. Clones a shape and translates one of its coordinates in a given direction. */
void Packing::expandShapeOnBC(const RSAShape *shape, double translation, size_t translateCoordIdx) {
    RSAShape *shapeClone = shape->clone();
    RSAVector trans;
    trans[translateCoordIdx] = translation;
    shapeClone->translate(trans);
    this->addShape(shapeClone);
}
