//--------------------------------------------------------------------------------------------
// Class representing vectors and operations on them
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include "DynamicVector.h"

#include <stdexcept>



// Nothing special here, let the DynamicMatrix class do the magic
//--------------------------------------------------------------------------------------------
DynamicVector::DynamicVector()
{

}

DynamicVector::DynamicVector(const DynamicVector & other) : v(other.v)
{
    
}

DynamicVector::DynamicVector(const DynamicVector && other) : v(std::move(other.v))
{
    
}

DynamicVector & DynamicVector::operator=(const DynamicVector & other)
{
    this->v = other.v;
    return *this;
}

DynamicVector & DynamicVector::operator=(DynamicVector && other)
{
    this->v = std::move(other.v);
    return *this;
}

DynamicVector::DynamicVector(std::initializer_list<double> _coords) : v(DynamicMatrix(_coords.size(), 1, _coords))
{
    
}

DynamicVector::DynamicVector(vsize_t dim, double *_arr) : v(DynamicMatrix(dim, 1, _arr))
{

}

DynamicVector::DynamicVector(vsize_t dim, double _fill) : v(DynamicMatrix(dim, 1, _fill))
{

}

DynamicVector::DynamicVector(vsize_t dim) : v(DynamicMatrix(dim, 1, 0.))
{

}

DynamicVector::DynamicVector(const DynamicMatrix & _v) : v(_v)
{

}

DynamicVector::DynamicVector(DynamicMatrix && _v) : v(std::move(_v))
{

}


// Operators.
//--------------------------------------------------------------------------------------------

DynamicVector operator+(const DynamicVector & _v1, const DynamicVector & _v2)
{
    return DynamicVector(std::move(_v1.v + _v2.v));
}

// Scalar product
//--------------------------------------------------------------------------------------------
double operator*(const DynamicVector & _v1, const DynamicVector & _v2)
{
    if (_v1.getDimension() != _v2.getDimension())
        throw std::invalid_argument("_v1.getDimension() != _v2.getDimension()");
    double prod = 0;
    for (DynamicVector::vsize_t i = 0; i < _v1.getDimension(); i++)
        prod += _v1(i) * _v2(i);
    return prod;
}

// Multiplication by scalar
//--------------------------------------------------------------------------------------------
DynamicVector operator*(double _x, DynamicVector _v)
{
    return DynamicVector(std::move(_x * _v.v));
}

// Multiplication by scalar (different operands order)
//--------------------------------------------------------------------------------------------
DynamicVector operator*(const DynamicVector & _v, double _x)
{
    return DynamicVector(std::move(_x * _v.v));
}

// Linear transformation (multiplication of matrix and this)
//--------------------------------------------------------------------------------------------
DynamicVector operator*(const DynamicMatrix & _m, const DynamicVector & _v)
{
    return DynamicVector(std::move(_m * _v.v));
}

DynamicVector operator-(const DynamicVector & _v1, const DynamicVector & _v2)
{
    return DynamicVector(std::move(_v1.v - _v2.v));
}

// Cross product
//--------------------------------------------------------------------------------------------
DynamicVector operator^(const DynamicVector & _v1, const DynamicVector & _v2)
{
    if (_v1.getDimension() != 3)
        throw std::runtime_error("_v1 not a R^3 vector");
    if (_v2.getDimension() != 3)
        throw std::runtime_error("_v2 not a R^3 vector");
    return DynamicVector({
        _v1(1) * _v2(2) - _v1(2) * _v2(1),
        _v1(2) * _v2(0) - _v1(0) * _v2(2),
        _v1(0) * _v2(1) - _v1(1) * _v2(0) });
}

bool operator==(const DynamicVector & _v1, const DynamicVector & _v2)
{
    return (_v1.v == _v2.v);
}

bool operator!=(const DynamicVector & _v1, const DynamicVector & _v2)
{
    return (_v1.v != _v2.v);
}

DynamicVector & DynamicVector::operator+=(const DynamicVector & other)
{
    this->v += other.v;
    return *this;
}

// Multiplication by scalar
//--------------------------------------------------------------------------------------------
DynamicVector & DynamicVector::operator*=(double x)
{
    this->v *= x;
    return *this;
}

DynamicVector & DynamicVector::operator-=(const DynamicVector & other)
{
    this->v -= other.v;
    return *this;
}

DynamicVector DynamicVector::operator-(void) const
{
    return DynamicVector(std::move(-this->v));
}

// Cross product assignment operator
//--------------------------------------------------------------------------------------------
DynamicVector & DynamicVector::operator^=(const DynamicVector & other)
{
    *this = std::move(*this ^ other);
    return *this;
}


double & DynamicVector::operator()(vsize_t coord)
{
    return this->v(coord, 0);
}

const double & DynamicVector::operator()(vsize_t coord) const
{
    return this->v(coord, 0);
}


// Print vector to _ostr output stream
//-----------------------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & _ostr, const DynamicVector & _v)
{
    _ostr << "(";
    for (DynamicVector::vsize_t i = 0; i < _v.v.getRows() - 1; i++)
        _ostr << _v.v(i, 0) << ", ";
    _ostr << _v.v(_v.v.getRows() - 1, 0) << ")";
    return _ostr;
}
