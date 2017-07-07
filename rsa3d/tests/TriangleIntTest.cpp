//--------------------------------------------------------------------------------------------
// Test of intersections::tri_tri3D algorithm. It generates output to be inspected manually
// in Mathematica
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//------------------------------------------------------------------------------------------

#include "../Intersection.h"
#include "../Vector.h"
#include "../RND.h"

#include <sstream>
#include <string>
#include <fstream>
#include <iostream>

namespace
{
    std::string triangleToWolfram(const Vector * _tri)
    {
        std::stringstream out;
        out << "Triangle[{{" << _tri[0](0) << ", " << _tri[0](1) << ", " << _tri[0](2) << "}, ";
        out << "{" << _tri[1](0) << ", " << _tri[1](1) << ", " << _tri[1](2) << "}, ";
        out << "{" << _tri[2](0) << ", " << _tri[2](1) << ", " << _tri[2](2) << "}}]";
        return out.str();
    }
}


void TriangleIntTest_run()
{
    RND rnd;
    std::ofstream out("triangles.nb", std::ofstream::out);
    Vector triangle1[] = {Vector(3), Vector(3), Vector(3)};
    Vector triangle2[] = {Vector(3), Vector(3), Vector(3)};
    int no = 0;
    
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                triangle1[j](k) = rnd.nextValue();
                triangle2[j](k) = rnd.nextValue();
            }
        }
        out << "tri1 = " << triangleToWolfram(triangle1) << ";" << std::endl;
        out << "tri2 = " << triangleToWolfram(triangle2) << ";" << std::endl;
        out << std::boolalpha;
        out << "Graphics3D[{tri1, tri2}, PlotLabel->\"intersection::tri_tri3D: " << intersection::tri_tri3D(triangle1, triangle2) << "; ";
        out << "pair no: " << no++ << "\"]" << std::endl;
        out << std::endl;
    }
    
    out.close();
    std::cout << ">> File triangles.nb generated" << std::endl;
}
