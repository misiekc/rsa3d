//--------------------------------------------------------------------------------------------
// Class representing positive oriented face of convex 3D polyhedron
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------


#ifndef _ORIENTED_FACE_H
    #define _ORIENTED_FACE_H


#include <vector>

#include "Vector.h"


class OrientedFace
{
private:
    typedef Vector<3> * vptr;
    
    std::vector<vptr> vertices;
    
public:
    OrientedFace(const std::vector<vptr> & _vertices);
    
    Vector<3>   getOrthogonal() const;
    Vector<3>   getNormal() const;
    double      pointSignedDistance(const Vector<3> & _point) const;
    double      pointSign(const Vector<3> & _point) const;
    void        translate(const Vector<3> & _translation);
    double      getSurface() const;
    bool        intersectsWith(const OrientedFace & face) const;
    
    // Returns vertices of this face (as vector of R^3 Vector)
    //----------------------------------------------------------------------------------------
    std::vector<vptr> getVertices() const
    {
        return this->vertices;
    }
};


#endif // _ORIENTED_FACE_H
