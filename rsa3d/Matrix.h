//-----------------------------------------------------------------------------------------------------------------------------
// Class representing matrices and operations on them
//-----------------------------------------------------------------------------------------------------------------------------
// (C)PKua 2017
//-----------------------------------------------------------------------------------------------------------------------------

#ifndef _MATRIX_H
    #define _MATRIX_H

#include <utility>
#include <ostream>
#include <string>

class Matrix
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
    Matrix();
    Matrix(const Matrix & other);
    Matrix(Matrix && other);
    Matrix & operator=(const Matrix & other);
    Matrix & operator=(Matrix && other);
    ~Matrix();
    
    // Other ctors
    Matrix(mxsize_t _rows, mxsize_t _cols, double _fill);
    Matrix(mxsize_t _rows, mxsize_t _cols, double **_arr);
    Matrix(mxsize_t _rows, mxsize_t _cols, double *_arr);
    Matrix(mxsize_t _rows, mxsize_t _cols, std::initializer_list<double> _arr);

    // Static functions that generate matrices
    static Matrix identity(mxsize_t _size);
    static Matrix rotation2D(double _a);
    static Matrix rotation3D(double _ax, double _ay, double _az);

    // Operators
    friend Matrix operator+(Matrix matrix1, const Matrix & matrix2);
    friend Matrix operator*(const Matrix & matrix1, const Matrix & matrix2);
    friend Matrix operator*(double x, Matrix matrix);
    friend Matrix operator*(Matrix matrix, double x);
    friend Matrix operator-(Matrix matrix1, const Matrix & matrix2);
    friend Matrix operator^(const Matrix & matrix1, const Matrix & matrix2);
    friend bool operator==(const Matrix & matrix1, const Matrix & matrix2);
    friend bool operator!=(const Matrix & matrix1, const Matrix & matrix2);

    Matrix & operator+=(const Matrix & other);
    Matrix & operator*=(const Matrix & other);
    Matrix & operator*=(double x);
    Matrix & operator-=(const Matrix & other);
    Matrix operator-(void) const;
    Matrix & operator^=(const Matrix & other);

    double & operator()(mxsize_t row, mxsize_t column);
    const double & operator()(mxsize_t row, mxsize_t column) const;
    friend std::ostream & operator<< (std::ostream & _stream, const Matrix & _matrix);

    // Other operations
    mxsize_t getRows() const;
    mxsize_t getCols() const;
    Matrix transpose() const;
    double det() const;
    Matrix inverse() const;
    double matrix_minor(mxsize_t _row, mxsize_t _column) const;
    std::string toString() const;
};

// Private inline function returning the element from the array of given address. Doesn't perform
// arguments validation
//-------------------------------------------------------------------------------------------------------
inline double & Matrix::_get(mxsize_t _row, mxsize_t _col)
{
    return arr[_row * cols + _col];
}

// Private inline function returning immutable element from the array of given address. Doesn't perform
// arguments validation
//-------------------------------------------------------------------------------------------------------
inline const double & Matrix::_get(mxsize_t _row, mxsize_t _col) const
{
    return arr[_row * cols + _col];
}

// Inline function returning number of rows in matrix
//-------------------------------------------------------------------------------------------------------
inline Matrix::mxsize_t Matrix::getRows() const
{
    return rows;
}

// Inline function returning number of columns in matrix
//-------------------------------------------------------------------------------------------------------
inline Matrix::mxsize_t Matrix::getCols() const
{
    return cols;
}

#endif  // _MATRIX_H
