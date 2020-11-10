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
    SurfaceFunction::MinMax valueSpan;

public:
    CurvedSurface(int dim, double s, double ndx, double vdx, std::unique_ptr<SurfaceFunction> surfaceFunction,
                  std::unique_ptr<RSABoundaryConditions> bc);

    /**
     * @brief See SurfaceFunction::calculateValueRange
     */
    [[nodiscard]] SurfaceFunction::MinMax calculateValueRange(const RSAVector &voxelPosition,
                                                              double voxelSpatialSize) const;

    /**
     * @brief SurfaceFunction::fillInLastCoordinate
     */
    void fillInLastCoordinate(RSAVector &position) const;

    [[nodiscard]] double getArea() const override;
    [[nodiscard]] SurfaceFunction::MinMax getValueSpan() const;
};


#endif //RSA3D_CURVEDSURFACE_H
