//--------------------------------------------------------------------------------------------
// Test of cuboid intersection algorithm. It uses standard triangle-triangle intersection
// algorithm to check results given by Cuboid::overlap
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------


#include <iostream>
#include <cassert>
#include <cstring>
#include <sstream>
#include <vector>

#include "../Vector.h"
#include "../Intersection.h"
#include "../shapes/Cuboid.h"
#include "../BoundaryConditions.h"
#include "../ShapeFactory.h"


// Helper methods
//--------------------------------------------------------------------------------------------
namespace
{
    Vector cuboid1_tris[12][3];
    Vector cuboid2_tris[12][3];
    const double sizex = 1;
    const double sizey = 1;
    const double sizez = 1;
    
    void obtain_tris(Cuboid * cuboid, Vector (&arr)[12][3]);
    
    
    // Dummy BoundaryConditions
    class MockBC : public BoundaryConditions
    {
        double distance2(double *p1, double *p2)
        {
            return 0;
        }
        
	    double * getTranslation(double *result, double* p1, double* p2)
	    {
	        result[0] = 0;
	        result[1] = 0;
	        result[2] = 0;
	        return result;
	    }
    };
} 


// Runs intersection::tri_tri3D selftest. Two simple almoast-intersecting and 
// almoast-not-intersecting cases are checked
//--------------------------------------------------------------------------------------------
void TriTriInt_selftest_run()
{
    Vector triangle1[] = {
        Vector{0, -1, 0},
        Vector{0, 1, 0},
        Vector{0, 0, 1}
    };
    
    // Triangle on z = 1.00000001 plane containing x = 0, y = 0
    Vector triangle2[] = {
        Vector{1.4, 0.8, 1.00000001},
        Vector{-1.3, 0, 1.00000001},
        Vector{0, 4.5, 1.00000001}
    };
    
    // Triangle on z = 0.99999999 plane containing x = 0, y = 0
    Vector triangle3[] = {
        Vector{1.4, 0.8, 0.99999999},
        Vector{-1.3, 0, 0.99999999},
        Vector{0, 4.5, 0.99999999}
    };
    
    Vector triangle4[] = {
        Vector{-1, -1, 0},
        Vector{1, 0, 2.00000001},
        Vector{-1, 1, 0}
    };
    
    Vector triangle5[] = {
        Vector{-1, -1, 0},
        Vector{1, 0, 1.99999999},
        Vector{-1, 1, 0}
    };
    
    std::cout << std::boolalpha;
    std::cout << ">> Performing intersection::tri_tri3D quick selftest..." << std::endl;
    assert(intersection::tri_tri3D(triangle1, triangle2) == false);
    assert(intersection::tri_tri3D(triangle1, triangle3) == true);
    assert(intersection::tri_tri3D(triangle1, triangle4) == false);
    assert(intersection::tri_tri3D(triangle1, triangle5) == true);
    assert(intersection::tri_tri3D(triangle2, triangle1) == false);
    assert(intersection::tri_tri3D(triangle3, triangle1) == true);
    assert(intersection::tri_tri3D(triangle4, triangle1) == false);
    assert(intersection::tri_tri3D(triangle5, triangle1) == true);
    std::cout << "Passed." << std::endl;
}


// Performs Cuboid::overlap algorithm check. It generates some random pairs of cuboids
// and compares result given by Cuboid::overlap and intersection::tri_tri3D
//--------------------------------------------------------------------------------------------
void CuboidIntTest_run()
{
    // Test intersection::tri_tri3D
    TriTriInt_selftest_run();
    
    RND rnd;
    std::stringstream stream;
    stream << "3 " << sizex << " " <<  sizey << " " << sizez;
	ShapeFactory::initShapeClass("Cuboid", stream.str());
	MockBC bc;
}


// Helper method. Obtains and saves triangles from cuboid's faces
//--------------------------------------------------------------------------------------------
void obtain_tris(Cuboid * cuboid, Vector (&arr)[12][3])
{
    Vector pos(3, cuboid->getPosition()); 
    Matrix orientation = cuboid->getOrientation();   
    Vector vert[] = {
        pos + orientation * Vector{ sizex / 2,  sizey / 2,  sizez / 2},
        pos + orientation * Vector{-sizex / 2,  sizey / 2,  sizez / 2},
        pos + orientation * Vector{ sizex / 2, -sizey / 2,  sizez / 2},
        pos + orientation * Vector{ sizex / 2,  sizey / 2, -sizez / 2},
        pos + orientation * Vector{ sizex / 2, -sizey / 2, -sizez / 2},
        pos + orientation * Vector{-sizex / 2,  sizey / 2, -sizez / 2},
        pos + orientation * Vector{-sizex / 2, -sizey / 2,  sizez / 2},
        pos + orientation * Vector{-sizex / 2, -sizey / 2, -sizez / 2} 
    };
    
    arr[0][0] = vert[0];
    arr[0][1] = vert[1];
    arr[0][2] = vert[2];

    arr[1][0] = vert[2];
    arr[1][1] = vert[6];
    arr[1][2] = vert[1];

    arr[2][0] = vert[0];
    arr[2][1] = vert[2];
    arr[2][2] = vert[3];

    arr[3][0] = vert[3];
    arr[3][1] = vert[2];
    arr[3][2] = vert[4];

    arr[4][0] = vert[7];
    arr[4][1] = vert[2];
    arr[4][2] = vert[6];

    arr[5][0] = vert[7];
    arr[5][1] = vert[4];
    arr[5][2] = vert[2];

    arr[6][0] = vert[1];
    arr[6][1] = vert[0];
    arr[6][2] = vert[3];

    arr[7][0] = vert[5];
    arr[7][1] = vert[3];
    arr[7][2] = vert[1];

    arr[8][0] = vert[7];
    arr[8][1] = vert[5];
    arr[8][2] = vert[1];

    arr[9][0] = vert[7];
    arr[9][1] = vert[1];
    arr[9][2] = vert[6];

    arr[10][0] = vert[7];
    arr[10][1] = vert[3];
    arr[10][2] = vert[5];

    arr[11][0] = vert[7];
    arr[11][1] = vert[4];
    arr[11][2] = vert[3];
}


