//--------------------------------------------------------------------------------------------
// Class representing positive oriented face of convex 3D polyhedron
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------


#include <vector>

#include "Vector.h"


class OrientedFace
{
private:
    std::vector<Vector<3> *> vertices;
    
public:
    OrientedFace(const std::vector<Vector<3> *> & _vertices);
    OrientedFace(Vector<3> ** _vertices, std::size_t _num);
    
    Vector<3>   getOrthogonal() const;
    Vector<3>   getNormal() const;
    double      pointSignedDistance(const Vector<3> & _point) const;
    double      pointSign(const Vector<3> & _point) const;
    void        translate(const Vector<3> & _translation);
    std::vector<Vector<3> *> getVertices() const;
};
