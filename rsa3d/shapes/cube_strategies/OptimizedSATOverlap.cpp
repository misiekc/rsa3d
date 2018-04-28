//
// Created by PKua on 05.12.17.
//

#include "OptimizedSATOverlap.h"

using std::abs;


// Optimized SAT algorithm from:
// http://www.jkh.me/files/tutorials/Separating%20Axis%20Theorem%20for%20Oriented%20Bounding%20Boxes.pdf

int OptimizedSATOverlap::overlap(const Shape<3, 0> *first, const Shape<3, 0> *second) const {
    auto cube1 = dynamic_cast<const Cuboid*>(first);
    auto cube2 = dynamic_cast<const Cuboid*>(second);

    Vector<3> position1(cube1->getPosition());
    Vector<3> position2(cube2->getPosition());
    Matrix<3, 3> orientation1 = cube1->getOrientation();
    Matrix<3, 3> orientation2 = cube2->getOrientation();

    // Centroid difference
    Vector<3> T = position2 - position1;

    double size[3];
    Cuboid::getSize(size);

    // Cuboid halfsizes
    double W = size[COORD::X] / 2;
    double H = size[COORD::Y] / 2;
    double D = size[COORD::Z] / 2;

    // Unit normal vectors
    Vector<3> Ax = orientation1 * Vector<3>{{1, 0, 0}};
    Vector<3> Ay = orientation1 * Vector<3>{{0, 1, 0}};
    Vector<3> Az = orientation1 * Vector<3>{{0, 0, 1}};
    Vector<3> Bx = orientation2 * Vector<3>{{1, 0, 0}};
    Vector<3> By = orientation2 * Vector<3>{{0, 1, 0}};
    Vector<3> Bz = orientation2 * Vector<3>{{0, 0, 1}};
    
    double Rxx = Ax * Bx;
    double Rxy = Ax * By;
    double Rxz = Ax * Bz;
    double Ryx = Ay * Bx;
    double Ryy = Ay * By;
    double Ryz = Ay * Bz;
    double Rzx = Az * Bx;
    double Rzy = Az * By;
    double Rzz = Az * Bz;

    // TODO: overlap detection between cuboids of different sizes
    if (abs(T * Ax) > W + abs(W * Rxx) + abs(H * Rxy) + abs(D * Rxz))           // CASE 1
        return false;
    else if (abs(T * Ay) > H + abs(W * Ryx) + abs(H * Ryy) + abs(D * Ryz))      // CASE 2
        return false;
    else if (abs(T * Az) > D + abs(W * Rzx) + abs(H * Rzy) + abs(D * Rzz))      // CASE 3
        return false;
    else if (abs(T * Bx) > abs(W * Rxx) + abs(H * Ryx) + abs(D * Rzx) + W)      // CASE 4
        return false;
    else if (abs(T * By) > abs(W * Rxy) + abs(H * Ryy) + abs(D * Rzy) + H)      // CASE 5
        return false;
    else if (abs(T * Bz) > abs(W * Rxz) + abs(H * Ryz) + abs(D * Rzz) + D)      // CASE 6
        return false;
    else if (abs(T * Az * Ryx - T * Ay * Rzx) > abs(H * Rzx) + abs(D * Ryx) + abs(H * Rxz) + abs(D * Rxy))      // CASE 7
        return false;
    else if (abs(T * Az * Ryy - T * Ay * Rzy) > abs(H * Rzy) + abs(D * Ryy) + abs(W * Rxz) + abs(D * Rxx))      // CASE 8
        return false;
    else if (abs(T * Az * Ryz - T * Ay * Rzz) > abs(H * Rzz) + abs(D * Ryz) + abs(W * Rxy) + abs(H * Rxx))      // CASE 9
        return false;
    else if (abs(T * Ax * Rzx - T * Az * Rxx) > abs(W * Rzx) + abs(D * Rxx) + abs(H * Ryz) + abs(D * Ryy))      // CASE 10
        return false;
    else if (abs(T * Ax * Rzy - T * Az * Rxy) > abs(W * Rzy) + abs(D * Rxy) + abs(W * Ryz) + abs(D * Ryx))      // CASE 11
        return false;
    else if (abs(T * Ax * Rzz - T * Az * Rxz) > abs(W * Rzz) + abs(D * Rxz) + abs(W * Ryy) + abs(H * Ryx))      // CASE 12
        return false;
    else if (abs(T * Ay * Rxx - T * Ax * Ryx) > abs(W * Ryx) + abs(H * Rxx) + abs(H * Rzz) + abs(D * Rzy))      // CASE 13
        return false;
    else if (abs(T * Ay * Rxy - T * Ax * Ryy) > abs(W * Ryy) + abs(H * Rxy) + abs(W * Rzz) + abs(D * Rzx))      // CASE 14
        return false;
    else if (abs(T * Ay * Rxz - T * Ax * Ryz) > abs(W * Ryz) + abs(H * Rxz) + abs(W * Rzy) + abs(H * Rzx))      // CASE 15
        return false;

    return true;
}

std::string OptimizedSATOverlap::getName() const {
    return "OptimizedSATOverlap";
}

void OptimizedSATOverlap::runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) const {
    Matrix<3, 3> orientation1 = cube1->getOrientation();
    Matrix<3, 3> orientation2 = cube2->getOrientation();
    orientation1 * Vector<3>{{1, 0, 0}};
    orientation1 * Vector<3>{{0, 1, 0}};
    orientation1 * Vector<3>{{0, 0, 1}};
    orientation2 * Vector<3>{{1, 0, 0}};
    orientation2 * Vector<3>{{0, 1, 0}};
    orientation2 * Vector<3>{{0, 0, 1}};
}
