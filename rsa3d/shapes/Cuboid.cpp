//----------------------------------------------------------------------------
// Class representing cuboid, derived after Shape class.
//----------------------------------------------------------------------------
// (C)PKua 2017
//----------------------------------------------------------------------------

#include "Cuboid.h"
#include "cube_strategies/MineOverlap.h"
#include "cube_strategies/OptimizedSATOverlap.h"

#include <iterator>
#include <algorithm>


// Static atributes declaration
//----------------------------------------------------------------------------

double          Cuboid::size[3];
double          Cuboid::neighbourListCellSize;
double          Cuboid::voxelSize;
double          Cuboid::minDimension;
Vector<3>       Cuboid::relativeVertices[VERTEX::NUM_OF];


OverlapStrategy * Cuboid::defaultStrategy = new OptimizedSATOverlap;
OverlapStrategy * Cuboid::strategy = Cuboid::defaultStrategy;


// Default constructor creating new Cuboid in (0, 0, 0) with size set in
// Cuboid::initClass
//----------------------------------------------------------------------------
Cuboid::Cuboid(const Matrix<3, 3> & orientation) : Shape<3, 0>(), orientation(orientation)
{

}


// Static method for class initialization
//----------------------------------------------------------------------------
// args - string arguments. Format: <size x> <size y> <size z>
//----------------------------------------------------------------------------
void Cuboid::initClass(const std::string &args)
{
    std::stringstream args_stream (args);
    int dimension;
    args_stream >> dimension;   // fetch dimension for backward compatibility

    // Fetch and assert dimensions
    std::copy(std::istream_iterator<double>(args_stream), std::istream_iterator<double>(),
              size);
    if (size[0] <= 0.0 || size[1] <= 0.0 || size[2] <= 0.0)
        throw std::runtime_error("Wrong Cuboid dimensions: " + args);

    // renormailze sizes to obtain unit volume
    double volume = size[0] * size[1] * size[2];
    double factor = 1.0 / pow(volume, 1./3.);
    std::transform(size, size + 3, size, [factor](double d) { return d * factor; });

    // Neighbour list cell size - diagonal length. It satisfies conditions at the least favourable case - Cuboids with
    // centers near opposite edges of a cell and with diagonals lying on the line connecting their centers
    neighbourListCellSize = std::sqrt(size[0] * size[0] + size[1] * size[1] + size[2] * size[2]);

    // Voxel size. It satisfies conditions at the least favourable case - center of Cuboid lies in the corner of a voxel
    // and voxel's diagonal is parallel to the smallest Cuboid edge
    minDimension = *std::min_element(size, size + 3) / 2;
    voxelSize = minDimension / std::sqrt(3);

    calculateRelativeVerties();
}

void Cuboid::calculateRelativeVerties() {
    relativeVertices[VERTEX::PPP] = Vector<3>{{ size[COORD::X] / 2,  size[COORD::Y] / 2,  size[COORD::Z] / 2}};
    relativeVertices[VERTEX::NPP] = Vector<3>{{-size[COORD::X] / 2,  size[COORD::Y] / 2,  size[COORD::Z] / 2}};
    relativeVertices[VERTEX::PNP] = Vector<3>{{ size[COORD::X] / 2, -size[COORD::Y] / 2,  size[COORD::Z] / 2}};
    relativeVertices[VERTEX::PPN] = Vector<3>{{ size[COORD::X] / 2,  size[COORD::Y] / 2, -size[COORD::Z] / 2}};
    relativeVertices[VERTEX::PNN] = Vector<3>{{ size[COORD::X] / 2, -size[COORD::Y] / 2, -size[COORD::Z] / 2}};
    relativeVertices[VERTEX::NPN] = Vector<3>{{-size[COORD::X] / 2,  size[COORD::Y] / 2, -size[COORD::Z] / 2}};
    relativeVertices[VERTEX::NNP] = Vector<3>{{-size[COORD::X] / 2, -size[COORD::Y] / 2,  size[COORD::Z] / 2}};
    relativeVertices[VERTEX::NNN] = Vector<3>{{-size[COORD::X] / 2, -size[COORD::Y] / 2, -size[COORD::Z] / 2}};

}

// Method creating (dynamically alocated) cuboid with random orientation.
// Used by ShapeFactory for shape generating
//----------------------------------------------------------------------------
Shape<3, 0> * Cuboid::create3D(RND *rnd)
{
    Cuboid * cuboid = new Cuboid(Matrix<3, 3>::rotation(
        rnd->nextValue() * 2 * M_PI,
        std::asin(rnd->nextValue() * 2 - 1),
        rnd->nextValue() * 2 * M_PI));

    return cuboid;
}


// Method creating (dynamically alocated) cuboid with random orientation with
// 2 degrees of freedom
//----------------------------------------------------------------------------
Shape<3, 0> * Cuboid::create2D(RND *rnd)
{
    Cuboid * cuboid = new Cuboid(Matrix<3, 3>::rotation(
        0,
        std::asin(rnd->nextValue() * 2 - 1),
        rnd->nextValue() * 2 * M_PI));

    return cuboid;
}


// Sets which overlap algorithm to use
//----------------------------------------------------------------------------
void Cuboid::setOverlapStrategy(OverlapStrategy * _strategy)
{
    strategy = _strategy;
}


// Returns overlap algorithm used
//----------------------------------------------------------------------------
OverlapStrategy * Cuboid::getOverlapStrategy()
{
    return Cuboid::strategy;
}


// Restores default overlap algorithm
//----------------------------------------------------------------------------
void Cuboid::restoreDefaultStrategy() {
    setOverlapStrategy(defaultStrategy);
}


// Returns neighbour list cell size determined during class initialization
//----------------------------------------------------------------------------
double Cuboid::getNeighbourListCellSize()
{
    return Cuboid::neighbourListCellSize;
}


// Returns initial voxel size determined during class initialization
//----------------------------------------------------------------------------
double Cuboid::getVoxelSize()
{
    return Cuboid::voxelSize;
}


// Checks whether there is an overlap between this and *s using chosen
// stategy
//----------------------------------------------------------------------------
int Cuboid::overlap(BoundaryConditions *bc, Shape *s)
{
    return strategy->overlap(this, (Cuboid*)s, bc);
}


// Checks whether given vertex (in this coordinates) lies in Cuboid
//----------------------------------------------------------------------------
bool Cuboid::pointInsideCuboid(const Vector<3> &vertex)
{
    for (unsigned short i = 0; i < 3; i++)
        if (std::abs(vertex[i]) > this->size[i] / 2)
            return false;
    return true;
}


// Returns volume of the Cuboid determined during class initialization
//----------------------------------------------------------------------------
double Cuboid::getVolume()
{
    return 1.;
}

// Checks whether the point with coordinates da lies inside excluded volume.
// It is an interior of a set of point which distanse from
//----------------------------------------------------------------------------
int Cuboid::pointInside(BoundaryConditions *bc, double* pos, double *orientation, double orientationRange)
{
    Vector<3> cuboidTranslation(this->position);
    Vector<3> pointTranslation(pos);

    // Transform point coordinates to Cuboid coordinate system
    double trans_arr[3];
    cuboidTranslation += Vector<3>(bc->getTranslation(trans_arr, this->position, pos));
    pointTranslation = this->orientation.transpose() * (pointTranslation - cuboidTranslation);

    double abs_point_coords[3];
    for (unsigned short i = 0; i < 3; i++)
        abs_point_coords[i] = std::abs(pointTranslation[i]);

    // Map which coords lie in this and which lie in this + minDimension
    bool liesInSmaller[3];
    bool liesInBigger[3];
    for (unsigned short i = 0; i < 3; i++) {
        liesInSmaller[i] = (abs_point_coords[i] <= this->size[i] / 2);
        liesInBigger[i] = (abs_point_coords[i] <= this->size[i] / 2 + minDimension);
    }

    // Check optimistic cases - "lies in this" and "doesn't lie in this + minDimension"
    if (!liesInBigger[0] || !liesInBigger[1] || !liesInBigger[2])       return false;
    if (liesInSmaller[0] && liesInSmaller[1] && liesInSmaller[2])       return true;

    // Check "pushed rectangles"
    if ((liesInSmaller[0] && liesInSmaller[1]) || (liesInSmaller[1] && liesInSmaller[2]) || (liesInSmaller[2] && liesInSmaller[0]))
        return true;

    // Check cylinders on edges
    if (liesInSmaller[0] && std::pow(abs_point_coords[1] - this->size[1] / 2, 2) + std::pow(abs_point_coords[2] - this->size[2] / 2, 2) <= std::pow(minDimension, 2))
        return true;
    if (liesInSmaller[1] && std::pow(abs_point_coords[2] - this->size[2] / 2, 2) + std::pow(abs_point_coords[0] - this->size[0] / 2, 2) <= std::pow(minDimension, 2))
        return true;
    if (liesInSmaller[2] && std::pow(abs_point_coords[0] - this->size[0] / 2, 2) + std::pow(abs_point_coords[1] - this->size[1] / 2, 2) <= std::pow(minDimension, 2))
        return true;

    // Check spheres in vertices
    if (std::pow(abs_point_coords[0] - this->size[0] / 2, 2) + std::pow(abs_point_coords[1] - this->size[1] / 2, 2) +
        std::pow(abs_point_coords[2] - this->size[2] / 2, 2) <= std::pow(minDimension, 2))
        return true;

    return false;
}

// Returns Cuboid orientation
//----------------------------------------------------------------------------
Matrix<3, 3> Cuboid::getOrientation() const
{
    return this->orientation;
}

// Returns Cuboid size
//----------------------------------------------------------------------------
double * Cuboid::getSize(double * arr)
{
    std::copy(size, size + 3, arr);
    return arr;
}

std::string Cuboid::toPovray() const
{
	std::string s = "";
	s += "  box { < ";
	for(unsigned short i=0; i<3; i++){
		s += std::to_string(-this->size[i]/2);
		if (i<2)
			s+= ", ";
	}
	s += ">, <";
	for(unsigned short i=0; i<3; i++){
		s += std::to_string(this->size[i]/2);
		if (i<2)
			s+= ", ";
	}
	s += ">\n";
	s += "    matrix < \n    ";
	for (unsigned short i=0; i<3; i++){
		for (unsigned short j=0; j<3; j++){
			s += std::to_string(this->orientation(j, i));
			if (j<2)
				s+= ", ";
		}
		s+= ",\n    ";
	}
	for(unsigned short i=0; i<3; i++){
		s += std::to_string(this->position[i]);
		if (i<2)
			s+= ", ";
	}
	s += "\n    >\n";
	s += "    texture { pigment { color Red } }\n  }\n";
	return s;
}

// Returns Wolfram Language for drawing shape
//----------------------------------------------------------------------------
std::string Cuboid::toWolfram() const
{
    std::stringstream out;
    out << "cube" << this->no << " = " << std::endl;
    out << "    GeometricTransformation[" << std::endl;
    out << "        Cuboid[{" << (-this->size[0] / 2) << ", " << (-this->size[1] / 2) << ", " << (-this->size[2] / 2) << "}, ";
    out << "{" << (this->size[0] / 2) << ", " << (this->size[1] / 2) << ", " << (this->size[2] / 2) << "}]," << std::endl;
    out << "        AffineTransform[" << std::endl;
    out << "            {{{" << this->orientation(0, 0) << ", " << this->orientation(0, 1) << ", " << this->orientation(0, 2) << "}," << std::endl;
    out << "            {" << this->orientation(1, 0) << ", " << this->orientation(1, 1) << ", " << this->orientation(1, 2) << "}," << std::endl;
    out << "            {" << this->orientation(2, 0) << ", " << this->orientation(2, 1) << ", " << this->orientation(2, 2) << "}}," << std::endl;
    out << "            {" << this->position[0] << ", " << this->position[1] << ", " << this->position[2] << "}}]];";
    return out.str();
}

void Cuboid::store(std::ostream &f) const
{
	Shape<3, 0>::store(f);
	double d;
	for (unsigned short i=0; i<3; i++){
		for (unsigned short j=0; j<3; j++){
			d = this->orientation(i, j);
			f.write((char *)(&d), sizeof(double));
		}
	}
}

void Cuboid::restore(std::istream &f)
{
	Shape<3, 0>::restore(f);
	double d;
	for (unsigned short i=0; i<3; i++){
		for (unsigned short j=0; j<3; j++){
			f.read((char *)&d, sizeof(double));
			this->orientation(i, j) = d;
		}
	}
}

void Cuboid::obtainVertices(Vector<3> (&vertices)[VERTEX::NUM_OF], const Vector<3> &translation) {
    Vector<3> pos(this->position);
    for (std::size_t i = 0; i < VERTEX::NUM_OF; i++)
        vertices[i] = pos + translation + orientation * Cuboid::getRelativeVertex(i);
}

// Returns vertex position as if the cuboid centroid was in the origin of
// coordinate system and was aligned with all axes. Check Cuboid::initClass
// for vertices order
//----------------------------------------------------------------------------
const Vector<3>  Cuboid::getRelativeVertex(std::size_t index) {
    return relativeVertices[index];
}

