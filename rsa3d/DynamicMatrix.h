//-----------------------------------------------------------------------------------------------------------------------------
// Class representing matrices and operations on them
//-----------------------------------------------------------------------------------------------------------------------------
// (C)PKua 2017
//-----------------------------------------------------------------------------------------------------------------------------

#ifndef _GENERAL_MATRIX_H
    #define _GENERAL_MATRIX_H

#include <utility>
#include <ostream>
#include <string>

class DynamicMatrix
{
public:
    typedef unsigned char mxsize_t;     // type for storing matrix size and indices
    typedef unsigned short mxsizel_t;   // type for storing mxsize_t variables product

private:
    double *arr;
    mxsize_t rows;
    mxsize_t cols;

    double & _get(mxsize_t row, mxsize_t column);
    const double & _get(mxsize_t row, mxsize_t column) const;

public:
    // Default, copy and move ctor, copy and move assingment operation, dtor
    DynamicMatrix();
    DynamicMatrix(const DynamicMatrix & other);
    DynamicMatrix(DynamicMatrix && other);
    DynamicMatrix & operator=(const DynamicMatrix & other);
    DynamicMatrix & operator=(DynamicMatrix && other);
    ~DynamicMatrix();
    
    // Other ctors
    DynamicMatrix(mxsize_t _rows, mxsize_t _cols, double _fill);
    DynamicMatrix(mxsize_t _rows, mxsize_t _cols, double **_arr);
    DynamicMatrix(mxsize_t _rows, mxsize_t _cols, double *_arr);
    DynamicMatrix(mxsize_t _rows, mxsize_t _cols, std::initializer_list<double> _arr);

    // Static functions that generate matrices
    static DynamicMatrix identity(mxsize_t _size);
    static DynamicMatrix rotation2D(double _a);
    static DynamicMatrix rotation3D(double _ax, double _ay, double _az);

    // Operators
    friend DynamicMatrix operator+(DynamicMatrix matrix1, const DynamicMatrix & matrix2);
    friend DynamicMatrix operator*(const DynamicMatrix & matrix1, const DynamicMatrix & matrix2);
    friend DynamicMatrix operator*(double x, DynamicMatrix matrix);
    friend DynamicMatrix operator*(DynamicMatrix matrix, double x);
    friend DynamicMatrix operator-(DynamicMatrix matrix1, const DynamicMatrix & matrix2);
    friend DynamicMatrix operator^(const DynamicMatrix & matrix1, const DynamicMatrix & matrix2);
    friend bool operator==(const DynamicMatrix & matrix1, const DynamicMatrix & matrix2);
    friend bool operator!=(const DynamicMatrix & matrix1, const DynamicMatrix & matrix2);

    DynamicMatrix & operator+=(const DynamicMatrix & other);
    DynamicMatrix & operator*=(const DynamicMatrix & other);
    DynamicMatrix & operator*=(double x);
    DynamicMatrix & operator-=(const DynamicMatrix & other);
    DynamicMatrix operator-(void) const;
    DynamicMatrix & operator^=(const DynamicMatrix & other);

    double & operator()(mxsize_t row, mxsize_t column);
    const double & operator()(mxsize_t row, mxsize_t column) const;
    friend std::ostream & operator<< (std::ostream & _stream, const DynamicMatrix & _matrix);

    // Other operations
    mxsize_t getRows() const;
    mxsize_t getCols() const;
    DynamicMatrix transpose() const;
    double det() const;
    DynamicMatrix inverse() const;
    double matrix_minor(mxsize_t _row, mxsize_t _column) const;
    std::string toString() const;
};

// Private inline function returning the element from the array of given address. Doesn't perform
// arguments validation
//-------------------------------------------------------------------------------------------------------
inline double & DynamicMatrix::_get(mxsize_t _row, mxsize_t _col)
{
    return arr[_row * cols + _col];
}

// Private inline function returning immutable element from the array of given address. Doesn't perform
// arguments validation
//-------------------------------------------------------------------------------------------------------
inline const double & DynamicMatrix::_get(mxsize_t _row, mxsize_t _col) const
{
    return arr[_row * cols + _col];
}

// Inline function returning number of rows in matrix
//-------------------------------------------------------------------------------------------------------
inline DynamicMatrix::mxsize_t DynamicMatrix::getRows() const
{
    return rows;
}

// Inline function returning number of columns in matrix
//-------------------------------------------------------------------------------------------------------
inline DynamicMatrix::mxsize_t DynamicMatrix::getCols() const
{
    return cols;
}

#endif  // _MATRIX_H
