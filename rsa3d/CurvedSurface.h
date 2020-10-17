//
// Created by Piotr Kubala on 17/10/2020.
//

#ifndef RSA3D_CURVEDSURFACE_H
#define RSA3D_CURVEDSURFACE_H


#include "geometry/Vector.h"
#include "Surface.h"

class CurvedSurface : public Surface {
public:
    struct MinMax {
        double min{};
        double max{};
    };

    MinMax calculateValueRange(const RSAVector &voxelPosition, double voxelSpatialSize);
    void fillInLastCoordinate(RSAVector &position);
};


#endif //RSA3D_CURVEDSURFACE_H
