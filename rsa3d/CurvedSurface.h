//
// Created by Piotr Kubala on 17/10/2020.
//

#ifndef RSA3D_CURVEDSURFACE_H
#define RSA3D_CURVEDSURFACE_H

#include <memory>

#include "geometry/Vector.h"
#include "Surface.h"
#include "SurfaceFunction.h"

/**
 * @brief A curved surface dictated by underlying SurfaceFunction
 */
class CurvedSurface : public Surface {
private:
    std::unique_ptr<SurfaceFunction> surfaceFunction;

public:
    CurvedSurface(int dim, double s, double ndx, double vdx, std::unique_ptr<SurfaceFunction> surfaceFunction,
                  std::unique_ptr<RSABoundaryConditions> bc)
            : Surface(dim, s, ndx, vdx, std::move(bc)), surfaceFunction(std::move(surfaceFunction))
    { }

    /**
     * @brief See SurfaceFunction::calculateValueRange
     */
    [[nodiscard]] SurfaceFunction::MinMax calculateValueRange(const RSAVector &voxelPosition,
                                                              double voxelSpatialSize) const;

    /**
     * @brief SurfaceFunction::fillInLastCoordinate
     */
    void fillInLastCoordinate(RSAVector &position) const;

    double getArea() const override;
};


#endif //RSA3D_CURVEDSURFACE_H
