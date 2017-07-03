//-----------------------------------------------------------------------------------------------------------------------------
// Class representing matrices and operations on them
//-----------------------------------------------------------------------------------------------------------------------------
// (C)PKua 2017
//-----------------------------------------------------------------------------------------------------------------------------

#include "Matrix.h"

//#define MATRIX_DEBUG

#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <cmath>

#ifdef MATRIX_DEBUG
    #include <iostream>
#endif



// Default ctor generating zero 1x1 matrix
//--------------------------------------------------------------------------------------------
Matrix::Matrix() :
    rows(1),
    cols(1)
{
    arr = new double[1];
    arr[0] = 0;
}

// Copy ctor
//--------------------------------------------------------------------------------------------
Matrix::Matrix(const Matrix & other) :
    rows(other.rows),
    cols(other.cols)
{
#ifdef MATRIX_DEBUG
    std::cout << "in matrix copy ctor" << std::endl;
#endif

    // Alloc an array and copy data
    unsigned int max = other.rows * other.cols;
    arr = new double[max];
    std::copy(other.arr, other.arr + max, arr);
}

// Move ctor
//--------------------------------------------------------------------------------------------
Matrix::Matrix(Matrix && other) :
    arr(other.arr),
    rows(other.rows),
    cols(other.cols)
{
#ifdef MATRIX_DEBUG
    std::cout << "in matrix move ctor" << std::endl;
#endif

    other.arr = nullptr;
}

// Copy assingment operator. Returns this reference for operations such as
// (matrix = getmatrix()).getWidth()
//--------------------------------------------------------------------------------------------
Matrix & Matrix::operator=(const Matrix & other)
{
#ifdef MATRIX_DEBUG
    std::cout << "in matrix copy assignment" << std::endl;
#endif

    // Self assingment, skip
    if (this == &other)
        return *this;

    // Alloc a new array of different size, if needed
    mxsize_t max = other.rows * other.cols;
    if (max != rows * cols) {
        delete[] arr;
        arr = new double[max];
    }
    rows = other.rows;
    cols = other.cols;

    // Copy data
    std::copy(other.arr, other.arr + max, arr);

    return *this;
}

// Move assignment operator. Returns this reference for operations such as
// (matrix = getmatrix()).getWidth()
//--------------------------------------------------------------------------------------------
Matrix & Matrix::operator=(Matrix && other)
{
#ifdef MATRIX_DEBUG
    std::cout << "in matrix move assignment" << std::endl;
#endif

    delete[] arr;
    rows = other.rows;
    cols = other.cols;
    arr = other.arr;
    other.arr = nullptr;

    return *this;
}

// Dtor
//--------------------------------------------------------------------------------------------
Matrix::~Matrix()
{
#ifdef MATRIX_DEBUG
    std::cout << "in matrix dtor" << std::endl;
#endif

    delete[] arr;
}

// Ctor creating matrix of given size. Fills it with given value
//--------------------------------------------------------------------------------------------
// _rows, _cols - matrix dimensions
// _fill - value to be filled with
//--------------------------------------------------------------------------------------------
Matrix::Matrix(mxsize_t _rows, mxsize_t _cols, double _fill = 0) :
    rows(_rows),
    cols(_cols)
{
#ifdef MATRIX_DEBUG
    std::cout << "in fill ctor" << std::endl;
#endif

    if (_rows == 0 || _cols == 0)    throw std::invalid_argument("_rows == 0 || _cols == 0");

    // Alloc an array and fill
    unsigned int max = _rows * _cols;
    arr = new double[max];
    std::fill(arr, arr + max, _fill);
}

// Ctor creating matrix from given two-dimentional array, where the element from the i-th row
// and the j-th column is in arr[i-1][j-1]. Data from array are copied
//--------------------------------------------------------------------------------------------
// _rows, _cols - matrix dimensions
// _arr - array with matrix elements
//--------------------------------------------------------------------------------------------
Matrix::Matrix(mxsize_t _rows, mxsize_t _cols, double **_arr) :
    rows(_rows),
    cols(_cols)
{
#ifdef MATRIX_DEBUG
    std::cout << "in double** ctor" << std::endl;
#endif

    if (_rows == 0 || _cols == 0)    throw std::invalid_argument("_rows == 0 || _cols == 0");

    // Alloc an array and copy elements from given
    arr = new double[_rows * _cols];
    unsigned int arr_index = 0;
    for (mxsize_t i = 0; i < _rows; i++)
        for (mxsize_t j = 0; j < _cols; j++)
            arr[arr_index++] = _arr[i][j];
}

// Ctor creating matrix from given one-dimentional array. When desired matrix is
// / 1 2 3 |
// | 4 5 6 |
// \ 7 8 9 |
// then the array should be: {1, 2, 3, 4, 5, 6, 7, 8, 9}
//--------------------------------------------------------------------------------------------
// _rows, _cols - matrix dimensions
// _arr - array with matrix elements
//--------------------------------------------------------------------------------------------
Matrix::Matrix(mxsize_t _rows, mxsize_t _cols, double *_arr) :
    rows(_rows),
    cols(_cols)
{
#ifdef MATRIX_DEBUG
    std::cout << "in double* ctor" << std::endl;
#endif

    if (_rows == 0 || _cols == 0)    throw std::invalid_argument("_rows == 0 || _cols == 0");

    // Alloc an array and copy elements from given
    unsigned int max = _rows * _cols;
    arr = new double[max];
    std::copy(_arr, _arr + max, arr);
}

// Ctor creating matrix from initializer lsit. When desired matrix is
// / 1 2 3 |
// | 4 5 6 |
// \ 7 8 9 |
// then the list should be: {1, 2, 3, 4, 5, 6, 7, 8, 9}
//--------------------------------------------------------------------------------------------
// _rows, _cols - matrix dimensions
// _arr - initializer lsit with matrix elements
//--------------------------------------------------------------------------------------------
Matrix::Matrix(mxsize_t _rows, mxsize_t _cols, std::initializer_list<double> _arr) :
    rows(_rows),
    cols(_cols)
{
#ifdef MATRIX_DEBUG
    std::cout << "in il ctor" << std::endl;
#endif

    if (_rows == 0 || _cols == 0)            throw std::invalid_argument("_rows == 0 || _cols == 0");
    if ((mxsize_t)_arr.size() != _rows * _cols)    throw std::invalid_argument("initializer list size doesn't mach given dimensions");

    // Alloc an array and copy elements from list
    arr = new double[_arr.size()];
    std::copy(_arr.begin(), _arr.end(), arr);
}

// Generates identity matrix of given size
//--------------------------------------------------------------------------------------------
// _size - matrix size
//--------------------------------------------------------------------------------------------
Matrix Matrix::identity(short _size)
{
    Matrix matrix(_size, _size);
    for (mxsize_t i = 0; i < _size; i++)
        matrix._get(i, i) = 1;
    return matrix;
}

// Generates two-dimentional rotate matrix
//--------------------------------------------------------------------------------------------
// _a - rotate angle
//--------------------------------------------------------------------------------------------
Matrix Matrix::rotation2D(double _a)
{
    double sin_a = sin(_a);
    double cos_a = cos(_a);
    Matrix matrix(2, 2);
    matrix.arr[0] = cos_a;
    matrix.arr[1] = -sin_a;
    matrix.arr[2] = sin_a;
    matrix.arr[3] = cos_a;
    return matrix;
}

// Generates three-dimentional rotate matrix. The rotations are performed about X, Y and Z 
// axis in mentioned order
//--------------------------------------------------------------------------------------------
// _ax - counter-clockwise rotation angle about X axis
// _ay - counter-clockwise rotation angle about Y axis
// _az - counter-clockwise rotation angle about Z axis
//--------------------------------------------------------------------------------------------
Matrix Matrix::rotation3D(double _ax, double _ay, double _az)
{
    double sin_ax = sin(_ax);
    double sin_ay = sin(_ay);
    double sin_az = sin(_az);
    double cos_ax = cos(_ax);
    double cos_ay = cos(_ay);
    double cos_az = cos(_az);

    Matrix matrix(3, 3);
    matrix.arr[0] = cos_ay * cos_az;
    matrix.arr[1] = sin_ax * sin_ay * cos_az - cos_ax * sin_az;
    matrix.arr[2] = cos_ax * sin_ay * cos_az + sin_ax * sin_az;

    matrix.arr[3] = cos_ay * sin_az;
    matrix.arr[4] = sin_ax * sin_ay * sin_az + cos_ax * cos_az;
    matrix.arr[5] = cos_ax * sin_ay * sin_az - sin_ax * cos_az;
    
    matrix.arr[6] = -sin_ay;
    matrix.arr[7] = sin_ax * cos_ay;
    matrix.arr[8] = cos_ax * cos_ay;

    return matrix;
}

// Addition of 2 matrices
//--------------------------------------------------------------------------------------------
Matrix operator+(Matrix matrix1, const Matrix & matrix2)
{
    return (matrix1 += matrix2);
}

// Multiplication of 2 matrices
//--------------------------------------------------------------------------------------------
Matrix operator*(const Matrix & matrix1, const Matrix & matrix2)
{
    if (matrix1.cols != matrix2.rows)
        throw std::invalid_argument("cannot multiply matrices; matrix1.getCols() != matrix2.getRows()");

    Matrix::mxsize_t kmax = matrix1.cols;
    Matrix ret(matrix1.rows, matrix2.cols);
    for (Matrix::mxsize_t i = 0; i < matrix1.rows; i++)
        for (Matrix::mxsize_t j = 0; j < matrix2.cols; j++)
            for (Matrix::mxsize_t k = 0; k < kmax; k++)
                ret._get(i, j) += matrix1._get(i, k) * matrix2._get(k, j);
    return ret;
}

// Multiplying matrix by scalar
//--------------------------------------------------------------------------------------------
Matrix operator*(double x, Matrix matrix)
{
    return (matrix *= x);
}

Matrix operator*(Matrix matrix, double x)
{
    return (matrix *= x);
}

// Subtraction of 2 matrices
//--------------------------------------------------------------------------------------------
Matrix operator-(Matrix matrix1, const Matrix & matrix2)
{
    return (matrix1 -= matrix2);
}

// Cross product operator between two R3 vector
//--------------------------------------------------------------------------------------------
Matrix operator^(const Matrix & matrix1, const Matrix & matrix2)
{
    if (matrix1.getRows() != 3 || matrix1.getCols() != 1)
        throw std::runtime_error("matrix1 not in M(3, 1, R)");
    if (matrix2.getRows() != 3 || matrix2.getCols() != 1)
        throw std::runtime_error("matrix2 not in M(3, 1, R)");
    return Matrix(3, 1, {
        matrix1._get(1, 0) * matrix2._get(2, 0) - matrix1._get(2, 0) * matrix2._get(1, 0),
        matrix1._get(2, 0) * matrix2._get(0, 0) - matrix1._get(0, 0) * matrix2._get(2, 0),
        matrix1._get(0, 0) * matrix2._get(1, 0) - matrix1._get(1, 0) * matrix2._get(0, 0) });
}

// Equality operator of 2 matrices
//--------------------------------------------------------------------------------------------
bool operator==(const Matrix & matrix1, const Matrix & matrix2)
{
    if (matrix1.rows != matrix2.rows || matrix1.cols != matrix2.cols)
        throw std::invalid_argument("cannot compare matrices of different sizes");

    unsigned int max = matrix1.rows * matrix1.cols;
    for (unsigned int i = 0; i < max; i++)
        if (matrix1.arr[i] != matrix2.arr[i])
            return false;
    return true;
}

// Inequality operator of 2 matrices
//--------------------------------------------------------------------------------------------
bool operator!=(const Matrix & matrix1, const Matrix & matrix2)
{
    return !(matrix1 == matrix2);
}

// Addition assignment operator
//--------------------------------------------------------------------------------------------
Matrix & Matrix::operator+=(const Matrix & other)
{
    if (rows != other.rows || cols != other.cols)
        throw std::invalid_argument("cannot add matrices of different sizes");

    unsigned int max = rows * cols;
    for (unsigned int i = 0; i < max; i++)
        arr[i] += other.arr[i];
    return *this;
}

// Matrix multiplication assignment operator
//--------------------------------------------------------------------------------------------
Matrix & Matrix::operator*=(const Matrix & other)
{
    *this = std::move(*this * other);
    return *this;
}

// Matrix multiplication by scalar assignment operator
//--------------------------------------------------------------------------------------------
Matrix & Matrix::operator*=(double x)
{
    if (x == 0) {
        std::fill(arr, arr + rows * cols, 0);
    } else if (x != 1) {
        unsigned int max = rows * cols;
        for (unsigned int i = 0; i < max; i++)
            arr[i] *= x;
    }
    return *this;
}

// Subtraction assignment operator
//--------------------------------------------------------------------------------------------
Matrix & Matrix::operator-=(const Matrix & other)
{
    if (rows != other.rows || cols != other.cols)
        throw new std::invalid_argument("cannot subtract matrices of different sizes");

    unsigned int max = rows * cols;
    for (unsigned int i = 0; i < max; i++)
        arr[i] -= other.arr[i];
    return *this;
}

// Cross product assignment operator between two R3 vectors
//--------------------------------------------------------------------------------------------
Matrix & Matrix::operator^=(const Matrix & other)
{
    *this = std::move(*this ^ other);
    return *this;
}

// Unary minus operator
//--------------------------------------------------------------------------------------------
Matrix Matrix::operator-(void) const
{
    Matrix ret(rows, cols);
    unsigned int max = rows * cols;
    for (unsigned int i = 0; i < max; i++)
        ret.arr[i] = -arr[i];
    return ret;
}


// R/W access oparator. In order to get the element from i-th row and j-th column, one should
// use matrix(i-1, j-1). Check if indexes are in bounds and throws std::invalid_argument
// exception if needed
//--------------------------------------------------------------------------------------------
double & Matrix::operator()(mxsize_t _row, mxsize_t _column)
{
    if (_row >= rows)       throw std::invalid_argument("row index out of bounds");
    if (_column >= cols)    throw std::invalid_argument("column index out of bounds");

    return _get(_row, _column);
}

// Read-only access oparator. In order to get the element from i-th row and j-th column, one
// should use matrix(i-1, j-1). Check if indexes are in bounds and throws
// std::invalid_argument exception if needed
//--------------------------------------------------------------------------------------------
const double & Matrix::operator()(mxsize_t _row, mxsize_t _column) const
{
    if (_row >= rows)       throw std::invalid_argument("row index out of bounds");
    if (_column >= cols)    throw std::invalid_argument("column index out of bounds");

    return _get(_row, _column);
}

// Returns transposition of the matrix
//--------------------------------------------------------------------------------------------
Matrix Matrix::transpose() const
{
    Matrix matrix(cols, rows);
    for (mxsize_t i = 0; i < rows; i++)
        for (mxsize_t j = 0; j < cols; j++)
            matrix._get(j, i) = _get(i, j);

    return matrix;
}

// Returns determinant of the matrix
//--------------------------------------------------------------------------------------------
double Matrix::det() const
{
    if (rows != cols)  throw std::runtime_error("rows != cols");

    // Small matrices
    if (rows == 1)  return arr[0];
    if (rows == 2)  return arr[0] * arr[3] - arr[1] * arr[2];
    if (rows == 3)  return arr[0] * arr[4] * arr[8] + arr[1] * arr[5] * arr[6] + arr[2] * arr[3] * arr[7]
        - arr[2] * arr[4] * arr[6] - arr[1] * arr[3] * arr[8] - arr[0] * arr[5] * arr[7];

    // Bigger matrices - Laplace expansion
    double determinant = 0.0;
    bool even_perm = true;
    for (mxsize_t i = 0; i < rows; i++) {
        if (even_perm)  determinant += _get(0, i) * matrix_minor(0, i);
        else            determinant -= _get(0, i) * matrix_minor(0, i);
        even_perm = !even_perm;
    }
    return determinant;
}

// Returns inversion of the matrix. Throws an exception, when matrix isn't square or isn't
// invertible
//--------------------------------------------------------------------------------------------
Matrix Matrix::inverse() const
{
    if (rows != cols)    throw std::runtime_error("rows != cols");

    // Small matrices
    if (rows == 1) {
        if (arr[0] == 0)  throw std::runtime_error("cannot invert, det == 0");
        return Matrix(1, 1, 1 / arr[0]);
    } else if (rows == 2) {
        double determinant = det();
        Matrix out(2, 2);
        if (determinant == 0)  throw std::runtime_error("cannot invert, det == 0");
        
        out.arr[0] = arr[3] / determinant;
        out.arr[1] = -arr[1] / determinant;
        out.arr[2] = -arr[2] / determinant;
        out.arr[3] = arr[0] / determinant;
        return out;
    }
    // Bigger matrices
    else {
        double determinant = det();
        Matrix out(rows, rows);
        if (determinant == 0)  throw std::runtime_error("cannot invert, det == 0");

        for (mxsize_t i = 0; i < rows; i++)
            for (mxsize_t j = 0; j < cols; j++)
                out._get(i, j) = ((i + j) % 2 == 0) ? (matrix_minor(j, i) / determinant) : (-matrix_minor(j, i) / determinant);
        return out;
    }
}

// Calculates minor of matrix created by removing (_row + 1) row and (_column + 1) column 
//--------------------------------------------------------------------------------------------
double Matrix::matrix_minor(mxsize_t _row, mxsize_t _column) const
{
    if (rows != cols)    throw std::runtime_error("rows != cols");
    if (rows == 1)        throw std::runtime_error("zero-size minor");
    if (rows == 2)        return arr[3 ^ (_row << 1 | _column)];

    Matrix minor_mx(rows - 1, rows - 1);

    if (_row > 0) {
        // Copy left-top part (if exists)
        if (_column > 0)
            for (mxsize_t i = 0; i < _row; i++)
                for (mxsize_t j = 0; j < _column; j++)
                    minor_mx._get(i, j) = _get(i, j);

        // Copy right-top part (if exists)
        if (_column < cols - 1)
            for (mxsize_t i = 0; i < _row; i++)
                for (mxsize_t j = _column + 1; j < cols; j++)
                    minor_mx._get(i, j - 1) = _get(i, j);
    }

    if (_row < rows - 1) {
        // Copy left-bottom part (if exists)
        if (_column > 0)
            for (mxsize_t i = _row + 1; i < rows; i++)
                for (mxsize_t j = 0; j < _column; j++)
                    minor_mx._get(i - 1, j) = _get(i, j);

        // Copy right-bottom part (if exists)
        if (_column < cols - 1)
            for (mxsize_t i = _row + 1; i < rows; i++)
                for (mxsize_t j = _column + 1; j < cols; j++)
                    minor_mx._get(i - 1, j - 1) = _get(i, j);
    }

    return minor_mx.det();
}

// Depicts the matrix in the form of a string
//--------------------------------------------------------------------------------------------
std::string Matrix::toString() const
{
    std::stringstream sstream;
    sstream << *this;
    return sstream.str();
}

// Friend stream insertion operator
//--------------------------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & _stream, const Matrix & _matrix)
{
    for (Matrix::mxsize_t i = 0; i < _matrix.getRows(); i++) {
        for (Matrix::mxsize_t j = 0; j < _matrix.getCols(); j++) {
            _stream << _matrix._get(i, j) << " ";
        }
        _stream << std::endl;
    }

    return _stream;
}
