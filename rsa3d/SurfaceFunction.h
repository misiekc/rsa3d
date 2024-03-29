//
// Created by Piotr Kubala on 17/10/2020.
//

#ifndef RSA3D_SURFACEFUNCTION_H
#define RSA3D_SURFACEFUNCTION_H

#include "geometry/Vector.h"
#include "utils/Utils.h"

/**
 * @brief A function representing curved surface.
 */
class SurfaceFunction {
public:
    struct MinMax {
        double min{};
        double max{};

        [[nodiscard]] double getSpan() const { return this->max - this->min; }

        friend std::ostream &operator<<(std::ostream &out, const MinMax &minMax) {
            return out << "[" << minMax.min << ", " << minMax.max << "]";
        }
    };

    virtual ~SurfaceFunction() = default;

    /**
     * @brief Calculates the minimal and maximal value of a function in the region given by all but last coordinates
     * of a voxel determined by @a voxelPosition and @a voxelSpatialSize.
     */
    [[nodiscard]] virtual MinMax calculateValueRange(const RSAVector &voxelPosition, double voxelSpatialSize) const = 0;

    /**
     * @brief Calulates the value of the function from all but last coordinates in @a position an stores it in the
     * last coordinate
     */
    virtual void fillInLastCoordinate(RSAVector &position) const = 0;

    /**
     * @brief Return the real area of a surface given by this SurfaceFunction for the domain of linear size @a size
     */
    [[nodiscard]] virtual double getArea(double size) const = 0;
};

#endif //RSA3D_SURFACEFUNCTION_H
