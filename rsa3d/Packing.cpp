//
// Created by PKua on 17.06.18.
//

#include "Packing.h"
#include "shape/ShapeFactory.h"
#include "utils/Utils.h"
#include <fstream>

Packing::~Packing() {
    for (auto shape : this->packing)
        delete shape;
}

Packing::Packing(const Packing &other) {
    for (auto &shape : other.packing)
        this->addShape(shape->clone());
}

Packing &Packing::operator=(Packing other) {
    std::swap(this->packing, other.packing);
    return *this;
}

Packing::Packing(Packing &&other) noexcept : packing(std::move(other.packing)) {
    other.packing.clear();
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

void Packing::expandOnPBC(double linearSize) {
    Expects(linearSize > 0.0);
    double circusphereRadius = RSAShape::getCircumsphereRadius();

    for (std::size_t i = 0; i < RSA_SPATIAL_DIMENSION; i++) {
        std::size_t oldSize = this->size();
        for (std::size_t j = 0; j < oldSize; j++) {
            auto shape = (*this)[j];
            RSAVector position = shape->getPosition();
            if (position[i] < circusphereRadius)
                expandShapeOnBC(shape, linearSize, i);
            else if (position[i] > linearSize - circusphereRadius)
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

void Packing::removeShape(std::size_t index) {
    delete (*this)[index];  // free range check
    this->packing.erase(this->packing.begin() + index);
}

const RSAShape *Packing::operator[](std::size_t index) const {
    if (index >= this->size())
        throw std::runtime_error("index >= size");
    return this->packing[index];
}

void Packing::clear() {
    for (auto shape : this->packing)
        delete shape;
    this->packing.clear();
}

double Packing::getParticlesVolume(unsigned short dim) const {
    return std::accumulate(this->packing.begin(), this->packing.end(), 0.0, [dim](auto sum, auto shape) {
        return sum + shape->getVolume(dim);
    });
}

void Packing::reducePackingVolume(double newPackingVolume, unsigned short dim) {
    std::size_t numberOfParticlesToKeep{};
    double partialVolume{};

    while (partialVolume < newPackingVolume && numberOfParticlesToKeep < this->packing.size()) {
        partialVolume += this->packing[numberOfParticlesToKeep]->getVolume(dim);
        numberOfParticlesToKeep++;
    }

    if (numberOfParticlesToKeep == 0) {
        this->clear();
    } else {
        std::for_each(this->packing.begin() + numberOfParticlesToKeep, this->packing.end(), [](auto shapePtr) {
            delete shapePtr;
        });
        this->packing.erase(this->packing.begin() + numberOfParticlesToKeep, this->packing.end());
    }
}
