//
// Created by PKua on 03.12.17.
//

#include "MineOverlap.h"

// Define this to disable edge-intersection fix in Cuboid::OverlapStrategy::MINE
//#define DISABLE_OVERLAP_FIX


namespace
{
    // Checks whether given segment intersects with axis-oriented rectangular face
    // with middle on the perpendicular axis (OP)
    //----------------------------------------------------------------------------
    // plane_pos - face-determined plane coordinate on OP
    // bound1 - positive face edge coordinate on the 1st axis perpendicular to OP (OB1)
    // bound2 - positive face edge coordinate on the 2nd axis perpendicular to OP (OB2)
    // p1p, p1b1, p1b2 - 1st segment point coordinates respectively on OP, OB1, OB2
    // p2p, p2b1, p2b2 - 2nd segment point coordinates respectively on OP, OB1, OB2
    //----------------------------------------------------------------------------
    inline bool checkSegmentFace(double plane_pos, double bound1, double bound2,
                                 double p1p, double p1b1, double p1b2,
                                 double p2p, double p2b1, double p2b2)
    {
        // Check weather a plane determined by face lies between segment points
        if ((plane_pos > p1p && plane_pos > p2p) || (plane_pos < p1p && plane_pos < p2p))
            return false;

        // Intersect plane's face with segment's line
        double b1i = ((p1b1 - p2b1) * plane_pos + p1p * p2b1 - p2p * p1b1) / (p1p - p2p);
        double b2i = ((p1b2 - p2b2) * plane_pos + p1p * p2b2 - p2p * p1b2) / (p1p - p2p);

        // Check whether intersection point lies on face
        if (std::abs(b1i) > bound1 || std::abs(b2i) > bound2)
            return false;
        else
            return true;
    }
}


bool MineOverlap::overlap(Cuboid *cube1, Cuboid *cube2, BoundaryConditions *bc) {
    // Prepare matrices of translations for operations on shapes;
    Vector<3> thisTranslation(cube1->getPosition());
    Vector<3> sTranslation(cube2->getPosition());
    double transArray[3];
    sTranslation += Vector<3>(bc->getTranslation(transArray, cube1->getPosition(), cube2->getPosition()));
    Matrix<3, 3> backwards_rot = cube1->getOrientation().transpose();

    // Transform s coordinates to this coordinate system
    Matrix<3, 3> new_orientation = backwards_rot * cube2->getOrientation();
    sTranslation = backwards_rot * (sTranslation - thisTranslation);

    Vector<3> v_trans[VERTEX::NUM_OF];    // Calculated vertices coordinates
    Vector<3> v_trans_bis[VERTEX::NUM_OF];    // Calculated vertices coordinates in swapped cuboid order

    double size[3];
    cube1->getSize(size);

    // Check whether vertices of s lie in this. TO OPTIMIZE
    for (size_t i = 0; i < VERTEX::NUM_OF; i++) {
        v_trans[i] = new_orientation * Cuboid::getRelativeVertex(i) + sTranslation;
        if (cube1->pointInsideCuboid(v_trans[i]))
            return true;
    }

    // Transform this coordinates to s coordinate system
    new_orientation = new_orientation.transpose();
    sTranslation = -(new_orientation * sTranslation);

    // Check whether vertices of this lie in s
    for (size_t i = 0; i < VERTEX::NUM_OF; i++) {
        v_trans_bis[i] = new_orientation * Cuboid::getRelativeVertex(i) + sTranslation;
        if (cube2->pointInsideCuboid(v_trans_bis[i]))
            return true;
    }

    // Check whether edges of s lie in this. TO OPTIMIZE
    if (checkSegment(cube1, v_trans[VERTEX::PPP], v_trans[VERTEX::PPN]) ||
        checkSegment(cube1, v_trans[VERTEX::PPN], v_trans[VERTEX::PNN]) ||
        checkSegment(cube1, v_trans[VERTEX::PNN], v_trans[VERTEX::PNP]) ||
        checkSegment(cube1, v_trans[VERTEX::PNP], v_trans[VERTEX::PPP]) ||
        checkSegment(cube1, v_trans[VERTEX::NNN], v_trans[VERTEX::NNP]) ||
        checkSegment(cube1, v_trans[VERTEX::NNP], v_trans[VERTEX::NPP]) ||
        checkSegment(cube1, v_trans[VERTEX::NPP], v_trans[VERTEX::NPN]) ||
        checkSegment(cube1, v_trans[VERTEX::NPN], v_trans[VERTEX::NNN]) ||
        checkSegment(cube1, v_trans[VERTEX::PPP], v_trans[VERTEX::NPP]) ||
        checkSegment(cube1, v_trans[VERTEX::PPN], v_trans[VERTEX::NPN]) ||
        checkSegment(cube1, v_trans[VERTEX::PNN], v_trans[VERTEX::NNN]) ||
        checkSegment(cube1, v_trans[VERTEX::PNP], v_trans[VERTEX::NNP]))
    {
        return true;
    }

#ifndef DISABLE_OVERLAP_FIX
    // Check whether edges of this lie in s. TO OPTIMIZE
    if (checkSegment(cube2, v_trans_bis[VERTEX::PPP], v_trans_bis[VERTEX::PPN]) ||
        checkSegment(cube2, v_trans_bis[VERTEX::PPN], v_trans_bis[VERTEX::PNN]) ||
        checkSegment(cube2, v_trans_bis[VERTEX::PNN], v_trans_bis[VERTEX::PNP]) ||
        checkSegment(cube2, v_trans_bis[VERTEX::PNP], v_trans_bis[VERTEX::PPP]) ||
        checkSegment(cube2, v_trans_bis[VERTEX::NNN], v_trans_bis[VERTEX::NNP]) ||
        checkSegment(cube2, v_trans_bis[VERTEX::NNP], v_trans_bis[VERTEX::NPP]) ||
        checkSegment(cube2, v_trans_bis[VERTEX::NPP], v_trans_bis[VERTEX::NPN]) ||
        checkSegment(cube2, v_trans_bis[VERTEX::NPN], v_trans_bis[VERTEX::NNN]) ||
        checkSegment(cube2, v_trans_bis[VERTEX::PPP], v_trans_bis[VERTEX::NPP]) ||
        checkSegment(cube2, v_trans_bis[VERTEX::PPN], v_trans_bis[VERTEX::NPN]) ||
        checkSegment(cube2, v_trans_bis[VERTEX::PNN], v_trans_bis[VERTEX::NNN]) ||
        checkSegment(cube2, v_trans_bis[VERTEX::PNP], v_trans_bis[VERTEX::NNP]))
    {
        return true;
    }
#endif

    return false;
}

// Checks whether a segment determined by point1 and point2 intersects with
// Cuboid
//----------------------------------------------------------------------------
bool MineOverlap::checkSegment(Cuboid *cube, const Vector<3> & point1, const Vector<3> & point2)
{
    double size[3];
    cube->getSize(size);
    double hsize_x = size[COORD::X] / 2;
    double hsize_y = size[COORD::Y] / 2;
    double hsize_z = size[COORD::Z] / 2;

    // Check intersections with all faces
    return checkSegmentFace(hsize_x, hsize_y, hsize_z, point1[COORD::X], point1[COORD::Y], point1[COORD::Z], point2[COORD::X], point2[COORD::Y], point2[COORD::Z]) ||
           checkSegmentFace(-hsize_x, hsize_y, hsize_z, point1[COORD::X], point1[COORD::Y], point1[COORD::Z], point2[COORD::X], point2[COORD::Y], point2[COORD::Z]) ||
           checkSegmentFace(hsize_y, hsize_z, hsize_x, point1[COORD::Y], point1[COORD::Z], point1[COORD::X], point2[COORD::Y], point2[COORD::Z], point2[COORD::X]) ||
           checkSegmentFace(-hsize_y, hsize_z, hsize_x, point1[COORD::Y], point1[COORD::Z], point1[COORD::X], point2[COORD::Y], point2[COORD::Z], point2[COORD::X]) ||
           checkSegmentFace(hsize_z, hsize_x, hsize_y, point1[COORD::Z], point1[COORD::X], point1[COORD::Y], point2[COORD::Z], point2[COORD::X], point2[COORD::Y]) ||
           checkSegmentFace(-hsize_z, hsize_x, hsize_y, point1[COORD::Z], point1[COORD::X], point1[COORD::Y], point2[COORD::Z], point2[COORD::X], point2[COORD::Y]);
}

std::string MineOverlap::getName() {
    return "MineOverlap";
}