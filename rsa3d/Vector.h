//--------------------------------------------------------------------------------------------
// Class representing vectors and operations on them
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _VECTOR_H
    #define _VECTOR_H

#include "Matrix.h"


class Vector
{
private:
    Matrix v;
    typedef Matrix::mxsize_t mxsize_t;
    
public:
    // Default, copy and move ctor, assignment operators
    Vector();
    Vector(const Vector & other);
    Vector(const Vector && other);
    Vector & operator=(const Vector & other);
    Vector & operator=(Vector && other);
    
    // Other ctors
    Vector(std::initializer_list<double> _coords);
    Vector(mxsize_t dim, double *_arr);
    Vector(mxsize_t dim, double _fill);
    explicit Vector(const Matrix & _v);
    explicit Vector(Matrix && _v);
    
    // Operators
    friend Vector operator+(const Vector & _v1, const Vector & _v2);
    friend double operator*(const Vector & _v1, const Vector & _v2);
    friend Vector operator*(double _x, Vector _v);
    friend Vector operator*(const Vector & _v, double _x);
    friend Vector operator*(const Matrix & _m, const Vector & _v);
    friend Vector operator-(const Vector & _v1, const Vector & _v2);
    friend Vector operator^(const Vector & _v1, const Vector & _v2);
    friend bool operator==(const Vector & _v1, const Vector & _v2);
    friend bool operator!=(const Vector & _v1, const Vector & _v2);

    Vector & operator+=(const Vector & other);
    Vector & operator*=(double x);
    Vector & operator-=(const Vector & other);
    Vector operator-(void) const;
    Vector & operator^=(const Vector & other);

    double & operator()(mxsize_t coord);
    const double & operator()(mxsize_t coord) const; 
    friend std::ostream & operator<<(std::ostream & _ostr, const Vector & _v);
    
    // Other operations
    Matrix::mxsize_t getDimension() const;
};

inline Matrix::mxsize_t Vector::getDimension() const
{
    return this->v.getRows();
}

#endif  // _VECTOR_H
