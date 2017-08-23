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
#include "../shapes/ConvexPolyhedron.h"
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
	    
	    Cuboid::setOverlapStrategy(Cuboid::OverlapStrategy::MINE);
	    int1 = (bool)cube1->overlap(&bc, cube2);
	    
	    Cuboid::setOverlapStrategy(Cuboid::OverlapStrategy::TRI_TRI);
	    int2 = (bool)cube1->overlap(&bc, cube2);
	    
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


// Performs Cuboid::overlap algorithm check. It generates some random pairs of cuboids
// and compares result given by Cuboid::overlap and ConvexPolyhedron::overlap
//--------------------------------------------------------------------------------------------
void CuboidIntSATTest_run()
{
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
	int     no = 0;
	
	std::cout << ">> Performing Cuboid::overlap and ConvexPolyhedron::overlap comparison..." << std::endl;
	for (int i = 0; i < tries; i++) {
	    cube1 = random_cuboid(&rnd);
	    cube2 = random_cuboid(&rnd);
	    
	    cube1->no = no;
	    no++;
	    
	    cube2->no = no;
	    no++;
	    
	    Cuboid::setOverlapStrategy(Cuboid::OverlapStrategy::MINE);
	    int1 = (bool)cube1->overlap(&bc, cube2);
	    
	    Cuboid::setOverlapStrategy(Cuboid::OverlapStrategy::SAT);
	    int2 = (bool)cube1->overlap(&bc, cube2);
	    
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
	nanoseconds     overlap_time, polyh_time, sat_time;
	
	std::cout << ">> Performing Cuboid::overlap strategy time comparison..." << std::endl;
	
	Cuboid::setOverlapStrategy(Cuboid::OverlapStrategy::MINE);
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
	
	Cuboid::setOverlapStrategy(Cuboid::OverlapStrategy::TRI_TRI);
	time_before = system_clock::now();
	for (int i = 0; i < tries; i++) {
	    cube1 = random_cuboid(&rnd);
	    cube2 = random_cuboid(&rnd);
	    
	    cube1->overlap(&bc, cube2);
	    
	    delete cube1;
	    delete cube2;
	}
	time_after = system_clock::now();
	polyh_time = duration_cast<std::chrono::nanoseconds>(time_after - time_before);
	polyh_time /= tries;

	Cuboid::setOverlapStrategy(Cuboid::OverlapStrategy::SAT);
	time_before = system_clock::now();
	for (int i = 0; i < tries; i++) {
	    cube1 = random_cuboid(&rnd);
	    cube2 = random_cuboid(&rnd);
	    
	    cube1->overlap(&bc, cube2);
	    
	    delete cube1;
	    delete cube2;
	}
	time_after = system_clock::now();
	sat_time = duration_cast<std::chrono::nanoseconds>(time_after - time_before);
	sat_time /= tries;
	
	std::cout << ">> Average custom algorithm time: " << overlap_time.count() << " ns" << std::endl;
	std::cout << ">> Average triangle algorithm time: " << polyh_time.count() << " ns" << std::endl;
	std::cout << ">> Average SAT time: " << sat_time.count() << " ns" << std::endl;
}
