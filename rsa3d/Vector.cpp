//--------------------------------------------------------------------------------------------
// Class representing vectors and operations on them
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include "Vector.h"

#include <stdexcept>



// Nothing special here, let the Matrix class do the magic
//--------------------------------------------------------------------------------------------
Vector::Vector()
{

}

Vector::Vector(const Vector & other) : v(other.v)
{
    
}

Vector::Vector(const Vector && other) : v(std::move(other.v))
{
    
}

Vector & Vector::operator=(const Vector & other)
{
    this->v = other.v;
    return *this;
}

Vector & Vector::operator=(Vector && other)
{
    this->v = std::move(other.v);
    return *this;
}

Vector::Vector(std::initializer_list<double> _coords) : v(Matrix(_coords.size(), 1, _coords))
{
    
}

Vector::Vector(mxsize_t dim, double *_arr) : v(Matrix(dim, 1, _arr))
{

}

Vector::Vector(mxsize_t dim, double _fill) : v(Matrix(dim, 1, _fill))
{

}

Vector::Vector(mxsize_t dim) : v(Matrix(dim, 1, 0.))
{

}

Vector::Vector(const Matrix & _v) : v(_v)
{

}

Vector::Vector(Matrix && _v) : v(std::move(_v))
{

}


// Operators.
//--------------------------------------------------------------------------------------------

Vector operator+(const Vector & _v1, const Vector & _v2)
{
    return Vector(std::move(_v1.v + _v2.v));
}

// Scalar product
//--------------------------------------------------------------------------------------------
double operator*(const Vector & _v1, const Vector & _v2)
{
    if (_v1.getDimension() != _v2.getDimension())
        throw std::invalid_argument("_v1.getDimension() != _v2.getDimension()");
    double prod = 0;
    for (Matrix::mxsize_t i = 0; i < _v1.getDimension(); i++)
        prod += _v1(i) * _v2(i);
    return prod;
}

// Multiplication by scalar
//--------------------------------------------------------------------------------------------
Vector operator*(double _x, Vector _v)
{
    return Vector(std::move(_x * _v.v));
}

// Multiplication by scalar (different operands order)
//--------------------------------------------------------------------------------------------
Vector operator*(const Vector & _v, double _x)
{
    return Vector(std::move(_x * _v.v));
}

// Linear transformation (multiplication of matrix and this)
//--------------------------------------------------------------------------------------------
Vector operator*(const Matrix & _m, const Vector & _v)
{
    return Vector(std::move(_m * _v.v));
}

Vector operator-(const Vector & _v1, const Vector & _v2)
{
    return Vector(std::move(_v1.v - _v2.v));
}

// Cross product
//--------------------------------------------------------------------------------------------
Vector operator^(const Vector & _v1, const Vector & _v2)
{
    if (_v1.getDimension() != 3)
        throw std::runtime_error("_v1 not a R^3 vector");
    if (_v2.getDimension() != 3)
        throw std::runtime_error("_v2 not a R^3 vector");
    return Vector({
        _v1(1) * _v2(2) - _v1(2) * _v2(1),
        _v1(2) * _v2(0) - _v1(0) * _v2(2),
        _v1(0) * _v2(1) - _v1(1) * _v2(0) });
}

bool operator==(const Vector & _v1, const Vector & _v2)
{
    return (_v1.v == _v2.v);
}

bool operator!=(const Vector & _v1, const Vector & _v2)
{
    return (_v1.v != _v2.v);
}

Vector & Vector::operator+=(const Vector & other)
{
    this->v += other.v;
    return *this;
}

// Multiplication by scalar
//--------------------------------------------------------------------------------------------
Vector & Vector::operator*=(double x)
{
    this->v *= x;
    return *this;
}

Vector & Vector::operator-=(const Vector & other)
{
    this->v -= other.v;
    return *this;
}

Vector Vector::operator-(void) const
{
    return Vector(std::move(-this->v));
}

// Cross product assignment operator
//--------------------------------------------------------------------------------------------
Vector & Vector::operator^=(const Vector & other)
{
    *this = std::move(*this ^ other);
    return *this;
}


double & Vector::operator()(mxsize_t coord)
{
    return this->v(coord, 0);
}

const double & Vector::operator()(mxsize_t coord) const
{
    return this->v(coord, 0);
}


// Print vector to _ostr output stream
//-----------------------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & _ostr, const Vector & _v)
{
    _ostr << "(";
    for (Vector::mxsize_t i = 0; i < _v.v.getRows() - 1; i++)
        _ostr << _v.v(i, 0) << ", ";
    _ostr << _v.v(_v.v.getRows() - 1, 0) << ")";
    return _ostr;
}
