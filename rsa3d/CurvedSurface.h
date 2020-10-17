//
// Created by Piotr Kubala on 17/10/2020.
//

#ifndef RSA3D_CURVEDSURFACE_H
#define RSA3D_CURVEDSURFACE_H

#include <memory>

#include "geometry/Vector.h"
#include "surfaces/NBoxPBC.h"
#include "SurfaceFunction.h"

class CurvedSurface : public NBoxPBC {
private:
    std::unique_ptr<SurfaceFunction> surfaceFunction;

public:
    CurvedSurface(int dim, double s, double ndx, double vdx, std::unique_ptr<SurfaceFunction> surfaceFunction)
            : NBoxPBC(dim, s, ndx, vdx), surfaceFunction(std::move(surfaceFunction))
    { }

    [[nodiscard]] SurfaceFunction::MinMax calculateValueRange(const RSAVector &voxelPosition,
                                                              double voxelSpatialSize) const;
    void fillInLastCoordinate(RSAVector &position) const;
};


#endif //RSA3D_CURVEDSURFACE_H
