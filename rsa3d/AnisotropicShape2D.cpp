//
// Created by PKua on 07.02.18.
//

#include "AnisotropicShape2D.h"

/**
 * Lazy pointInside - it uses angle-range point inside with full algle
 * @param bc BoundaryConditions used
 * @param da point to check
 * @return 0 if not inside
 */
int AnisotropicShape2D::pointInside(BoundaryConditions *bc, double *da) {
    return pointInside(bc, da, 0, 2 * M_PI);
}
