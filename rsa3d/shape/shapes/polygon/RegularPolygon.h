//
// Created by Michal Ciesla on 10.01.2023.
//

#ifndef RSA3D_REGULARPOLYGON_H
#define RSA3D_REGULARPOLYGON_H

#include "Polygon.h"

class RegularPolygon : public Polygon {
public:
    /**
     * @param args contains information about rounded polygon. See RegularDiskopolygonAttributes for format.
     */
    static void initClass(const std::string &attr);
};


#endif //RSA3D_REGULARPOLYGON_H
