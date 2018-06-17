//
// Created by PKua on 17.06.18.
//

#include "Packing.h"
#include "shape/ShapeFactory.h"
#include <fstream>

Packing::Packing(const Packing &other) {
    for (auto &shape : other.packing)
        this->addShape(shape->clone());
}

Packing &Packing::operator=(const Packing &other) {
    this->packing.clear();
    for (auto &shape : other.packing)
        this->addShape(shape->clone());
    return *this;
}

void Packing::addShape(const RSAShape *shape) {
    this->packing.push_back(shape_ptr(shape));
}

void Packing::removeShape(std::size_t index) {
    this->packing.erase(this->packing.begin() + index);
}

std::size_t Packing::getSize() const {
    return this->packing.size();
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
    this->removeShape(this->getSize() - 1);
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

const RSAShape *Packing::operator[](std::size_t index) const {
    return this->packing[index].get();
}

void Packing::expandOnPBC(double size, double expandMargin) {
    for (std::size_t i = 0; i < RSA_SPATIAL_DIMENSION; i++) {
        for (std::size_t j = 0; j < this->getSize(); j++) {
            auto shape = (*this)[j];
            const double *position = shape->getPosition();
            if (position[i] < expandMargin * size)
                expandShapeOnBC(shape, size, i);
            else if (position[i] > (1 - expandMargin) * size)
                expandShapeOnBC(shape, -size, i);
        }
    }
}

/* Helper method. Clones a shape and translates one of its coordinates in a given direction. */
void Packing::expandShapeOnBC(const RSAShape *shape, double translation, size_t translateCoordIdx) {
    RSAShape *shapeClone = shape->clone();
    std::array<double, RSA_SPATIAL_DIMENSION> trans{};
    trans.fill(0);
    trans[translateCoordIdx] = translation;
    shapeClone->translate(trans.data());
    this->addShape(shapeClone);
}
