//
// Created by PKua on 17.06.18.
//

#ifndef RSA3D_PACKING_H
#define RSA3D_PACKING_H

#include "shape/Shape.h"
#include <memory>

class Packing {
private:
    using shape_ptr = std::unique_ptr<const RSAShape>;

    std::vector<shape_ptr> packing;

    void expandShapeOnBC(const RSAShape *shape, double translation, size_t translateCoordIdx);

public:
    Packing(const Packing &other);
    Packing &operator=(const Packing &other);

    Packing() = default;
    virtual ~Packing() = default;
    Packing(Packing &&other) noexcept = default;
    Packing &operator=(Packing &&other) noexcept = default;

    void addShape(const RSAShape *shape);
    const RSAShape *operator[](std::size_t index) const;
    void removeShape(std::size_t index);
    std::size_t getSize() const;
    void store(std::ostream &out) const;
    void store(const std::string &filename) const;
    void restore(std::istream &in);
    void restore(const std::string &filename);

    /**
     * @brief Clones and translates shapes distant from the surface boundary not more than @a expandMargin x @a size
     * according to periodic boundary conditions
     * @param packing a packing to expand
     * @param size the size of the packing
     * @param expandMargin the distance (from 0 to 0.5) relative to @a size from the surface boundary from which shapes
     * will be translated
     */
    void expandOnPBC(double size, double expandMargin);
};


#endif //RSA3D_PACKING_H
