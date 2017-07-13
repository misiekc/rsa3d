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
#include <chrono>

#include "../Vector.h"
#include "../Intersection.h"
#include "../shapes/Cuboid.h"
#include "../BoundaryConditions.h"
#include "../ShapeFactory.h"

using namespace std::chrono;

// Helper methods
//--------------------------------------------------------------------------------------------
namespace
{
    const double    sizex = 1;
    const double    sizey = 1;
    const double    sizez = 1;
    const double    box_halfsize = 1;
    const int       tries = 10000;
    
	int             no = 0; 
    Vector<3>       cuboid1_tris[12][3];
    Vector<3>       cuboid2_tris[12][3];
    
    // Helper method. Obtains and saves triangles from cuboid's faces
    //--------------------------------------------------------------------------------------------
    void obtain_tris(Cuboid * cuboid, Vector<3> (&arr)[12][3])
    {
        Vector<3> pos(cuboid->getPosition()); 
        Matrix<3, 3> orientation = cuboid->getOrientation();   
        Vector<3> vert[] = {
            pos + orientation * Vector<3>{{ sizex / 2,  sizey / 2,  sizez / 2}},
            pos + orientation * Vector<3>{{-sizex / 2,  sizey / 2,  sizez / 2}},
            pos + orientation * Vector<3>{{ sizex / 2, -sizey / 2,  sizez / 2}},
            pos + orientation * Vector<3>{{ sizex / 2,  sizey / 2, -sizez / 2}},
            pos + orientation * Vector<3>{{ sizex / 2, -sizey / 2, -sizez / 2}},
            pos + orientation * Vector<3>{{-sizex / 2,  sizey / 2, -sizez / 2}},
            pos + orientation * Vector<3>{{-sizex / 2, -sizey / 2,  sizez / 2}},
            pos + orientation * Vector<3>{{-sizex / 2, -sizey / 2, -sizez / 2}} 
        };
        
        arr[0][0] = vert[0];
        arr[0][1] = vert[1];
        arr[0][2] = vert[2];

        arr[1][0] = vert[2];
        arr[1][1] = vert[1];
        arr[1][2] = vert[6];

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
        arr[7][1] = vert[1];
        arr[7][2] = vert[3];

        arr[8][0] = vert[7];
        arr[8][1] = vert[1];
        arr[8][2] = vert[5];

        arr[9][0] = vert[7];
        arr[9][1] = vert[6];
        arr[9][2] = vert[1];

        arr[10][0] = vert[7];
        arr[10][1] = vert[5];
        arr[10][2] = vert[3];

        arr[11][0] = vert[7];
        arr[11][1] = vert[3];
        arr[11][2] = vert[4];
    }
    
    // Helper method. Creates random Cuboid based on global parameters
    //--------------------------------------------------------------------------------------------
    Cuboid * random_cuboid(RND * rnd)
    {
	    double trans[3];
	    Cuboid * cube = (Cuboid*)ShapeFactory::createShape(rnd);
	    cube->no = no++;
	    trans[0] = (rnd->nextValue() * 2 - 1) * box_halfsize;
	    trans[1] = (rnd->nextValue() * 2 - 1) * box_halfsize;
	    trans[2] = (rnd->nextValue() * 2 - 1) * box_halfsize;
	    cube->translate(trans);
	    return cube;
    }
    
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
    Vector<3> triangle1[] = {
        Vector<3>{{0, -1, 0}},
        Vector<3>{{0, 1, 0}},
        Vector<3>{{0, 0, 1}}
    };
    
    // Triangle on z = 1.00000001 plane containing x = 0, y = 0
    Vector<3> triangle2[] = {
        Vector<3>{{1.4, 0.8, 1.00000001}},
        Vector<3>{{-1.3, 0, 1.00000001}},
        Vector<3>{{0, -4.5, 1.00000001}}
    };
    
    // Triangle on z = 0.99999999 plane containing x = 0, y = 0
    Vector<3> triangle3[] = {
        Vector<3>{{1.4, 0.8, 0.99999999}},
        Vector<3>{{-1.3, 0, 0.99999999}},
        Vector<3>{{0, -4.5, 0.99999999}}
    };
    
    Vector<3> triangle4[] = {
        Vector<3>{{-1, -1, 0}},
        Vector<3>{{1, 0, 2.00000001}},
        Vector<3>{{-1, 1, 0}}
    };
    
    Vector<3> triangle5[] = {
        Vector<3>{{-1, -1, 0}},
        Vector<3>{{1, 0, 1.99999999}},
        Vector<3>{{-1, 1, 0}}
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
	
	std::cout << std::boolalpha;
	Cuboid  *cube1, *cube2;
	bool    int1, int2;
	int     intersected = 0;
	int     disjunct = 0;
	int     missed = 0;
	
	std::cout << ">> Performing Cuboid::overlap and intersection::polyh_polyh comparison..." << std::endl;
	for (int i = 0; i < tries; i++) {
	    cube1 = random_cuboid(&rnd);
	    cube2 = random_cuboid(&rnd);
	    obtain_tris(cube1, cuboid1_tris);
	    obtain_tris(cube2, cuboid2_tris);
	    
	    int1 = (bool)cube1->overlap(&bc, cube2);
	    int2 = (bool)intersection::polyh_polyh(cuboid1_tris, 12, cuboid2_tris, 12);
	    if (int1 != int2) {
	        missed++;
	    } else {
	        if (int1)
	            intersected++;
	        else
	            disjunct++;
	    }
	    
	    if ((i % 10000) == 9999)
	        std::cout << (i + 1) << " pairs tested..." << std::endl;
	    
	    delete cube1;
	    delete cube2;
	}
	std::cout << ">> " << missed << " from " << tries << " intersection results missed" << std::endl;
	std::cout << ">> " << intersected << " cuboids overlapped, " << disjunct << " cuboids were disjunctive" << std::endl;
}

// Performs Cuboid::overlap and intersection::polyh_polyh time comparison
//--------------------------------------------------------------------------------------------
void CuboidIntTimeTest_run()
{
    RND rnd;
    std::stringstream stream;
    stream << "3 " << sizex << " " <<  sizey << " " << sizez;
	ShapeFactory::initShapeClass("Cuboid", stream.str());
	MockBC bc;
	
	std::cout << std::boolalpha;
	
	Cuboid  *cube1, *cube2;
	system_clock::time_point time_before, time_after;
	nanoseconds     overlap_time, polyh_time;
	
	std::cout << ">> Performing Cuboid::overlap and intersection::polyh_polyh time comparison..." << std::endl;
	time_before = system_clock::now();
	for (int i = 0; i < tries; i++) {
	    cube1 = random_cuboid(&rnd);
	    cube2 = random_cuboid(&rnd);
	    
	    cube1->overlap(&bc, cube2);
	    
	    delete cube1;
	    delete cube2;
	}
	time_after = system_clock::now();
	overlap_time = duration_cast<std::chrono::nanoseconds>(time_after - time_before);
	overlap_time /= tries;
	
	time_before = system_clock::now();
	for (int i = 0; i < tries; i++) {
	    cube1 = random_cuboid(&rnd);
	    cube2 = random_cuboid(&rnd);
	    obtain_tris(cube1, cuboid1_tris);
	    obtain_tris(cube2, cuboid2_tris);
	    
	    intersection::polyh_polyh(cuboid1_tris, 12, cuboid2_tris, 12);
	    
	    delete cube1;
	    delete cube2;
	}
	time_after = system_clock::now();
	polyh_time = duration_cast<std::chrono::nanoseconds>(time_after - time_before);
	polyh_time /= tries;
	
	std::cout << ">> Average Cuboid::overlap time: " << overlap_time.count() << " ns" << std::endl;
	std::cout << ">> Average intersection::polyh_polyh time: " << polyh_time.count() << " ns" << std::endl;
}
