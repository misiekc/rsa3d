//
// Created by PKua on 03.12.17.
//

#include "TriTriOverlap.h"
#include "../../Intersection.h"


namespace
{
    // Helper arrays
    Vector<3>       cuboid1_tris[12][3];
    Vector<3>       cuboid2_tris[12][3];
}

bool TriTriOverlap::overlap(const Cuboid *cube1, const Cuboid *cube2, BoundaryConditions *bc) {
    double trans_arr[3];
    Vector<3> translation(bc->getTranslation(trans_arr, cube1->getPosition(), cube2->getPosition()));

    obtainTris(cube1, cuboid1_tris, Vector<3>());
    obtainTris(cube2, cuboid2_tris, translation);

    return intersection::polyh_polyh(cuboid1_tris, 12, cuboid2_tris, 12);
}

// Helper method. Obtains and saves triangles from cuboid's faces
//--------------------------------------------------------------------------------------------
void TriTriOverlap::obtainTris(const Cuboid *cube, Vector<3> (&arr)[12][3], const Vector<3> &translation)
{
    double size[3];
    cube->getSize(size);
    Vector<3> vert[VERTEX::NUM_OF];
    cube->obtainVertices(vert, translation);

    arr[0][0] = vert[0];    arr[0][1] = vert[1];    arr[0][2] = vert[2];
    arr[1][0] = vert[2];    arr[1][1] = vert[1];    arr[1][2] = vert[6];
    arr[2][0] = vert[0];    arr[2][1] = vert[2];    arr[2][2] = vert[3];
    arr[3][0] = vert[3];    arr[3][1] = vert[2];    arr[3][2] = vert[4];
    arr[4][0] = vert[7];    arr[4][1] = vert[2];    arr[4][2] = vert[6];
    arr[5][0] = vert[7];    arr[5][1] = vert[4];    arr[5][2] = vert[2];
    arr[6][0] = vert[1];    arr[6][1] = vert[0];    arr[6][2] = vert[3];
    arr[7][0] = vert[5];    arr[7][1] = vert[1];    arr[7][2] = vert[3];
    arr[8][0] = vert[7];    arr[8][1] = vert[1];    arr[8][2] = vert[5];
    arr[9][0] = vert[7];    arr[9][1] = vert[6];    arr[9][2] = vert[1];
    arr[10][0] = vert[7];   arr[10][1] = vert[5];   arr[10][2] = vert[3];
    arr[11][0] = vert[7];   arr[11][1] = vert[3];   arr[11][2] = vert[4];
}

std::string TriTriOverlap::getName() {
    return "TriTriOverlap";
}

void TriTriOverlap::runOverheadOperations(const Cuboid *cube1, const Cuboid *cube2) {
    Vector<3> vert[VERTEX::NUM_OF];
    cube1->obtainVertices(vert, Vector<3>());
    cube2->obtainVertices(vert, Vector<3>());
}
