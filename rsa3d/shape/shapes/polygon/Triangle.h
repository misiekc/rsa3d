//
// Created by pkua on 10.04.19.
//

#ifndef RSA3D_TRIANGLE_H
#define RSA3D_TRIANGLE_H


#include "Polygon.h"

class Triangle : public Polygon {
public:
    ~Triangle() override = default;

    static void initClass(const std::string &args);

    /**
     * @brief Calculates the coordinates of the third triangle vertex assuming that its base is on X axis.
     *
     * It assumes that positions of two first vertices are (-c/2, 0), (c/2, 0) and the third one has y > 0
     * @param a the length of the right side
     * @param b the length of the left side
     * @param c the length of rhe horizontal side lying on x axis
     * @return the coordinates of third triangle vertex
     */
    static Vector<2> calculateThirdVertex(double a, double b, double c);
};


#endif //RSA3D_TRIANGLE_H
