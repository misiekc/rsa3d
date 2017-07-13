//--------------------------------------------------------------------------------------------
// Class representing vectors and operations on them
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _DYNAMIC_VECTOR_H
    #define _DYNAMIC_VECTOR_H

#include "DynamicMatrix.h"


class DynamicVector
{
private:
    DynamicMatrix v;
    typedef DynamicMatrix::mxsize_t vsize_t;
    
public:
    // Default, copy and move ctor, assignment operators
    DynamicVector();
    DynamicVector(const DynamicVector & other);
    DynamicVector(const DynamicVector && other);
    DynamicVector & operator=(const DynamicVector & other);
    DynamicVector & operator=(DynamicVector && other);
    
    // Other ctors
    DynamicVector(std::initializer_list<double> _coords);
    DynamicVector(vsize_t dim, double *_arr);
    DynamicVector(vsize_t dim, double _fill);
    DynamicVector(vsize_t dim);
    explicit DynamicVector(const DynamicMatrix & _v);
    explicit DynamicVector(DynamicMatrix && _v);
    
    // Operators
    friend DynamicVector operator+(const DynamicVector & _v1, const DynamicVector & _v2);
    friend double operator*(const DynamicVector & _v1, const DynamicVector & _v2);
    friend DynamicVector operator*(double _x, DynamicVector _v);
    friend DynamicVector operator*(const DynamicVector & _v, double _x);
    friend DynamicVector operator*(const DynamicMatrix & _m, const DynamicVector & _v);
    friend DynamicVector operator-(const DynamicVector & _v1, const DynamicVector & _v2);
    friend DynamicVector operator^(const DynamicVector & _v1, const DynamicVector & _v2);
    friend bool operator==(const DynamicVector & _v1, const DynamicVector & _v2);
    friend bool operator!=(const DynamicVector & _v1, const DynamicVector & _v2);

    DynamicVector & operator+=(const DynamicVector & other);
    DynamicVector & operator*=(double x);
    DynamicVector & operator-=(const DynamicVector & other);
    DynamicVector operator-(void) const;
    DynamicVector & operator^=(const DynamicVector & other);

    double & operator()(vsize_t coord);
    const double & operator()(vsize_t coord) const; 
    friend std::ostream & operator<<(std::ostream & _ostr, const DynamicVector & _v);
    
    // Other operations
    vsize_t getDimension() const;
};

inline DynamicVector::vsize_t DynamicVector::getDimension() const
{
    return this->v.getRows();
}

#endif  // _VECTOR_H
