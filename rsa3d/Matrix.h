//-----------------------------------------------------------------------------------------------------------------------------
// Class representing matrices and operations on them
//-----------------------------------------------------------------------------------------------------------------------------
// (C)PKua 2017
//-----------------------------------------------------------------------------------------------------------------------------

#ifndef _MATRIX_H
    #define _MATRIX_H

#include <array>
#include <utility>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>


// Forward declare the Vector class for friendship
//--------------------------------------------------------------------------------------------
template <std::size_t DIM, typename E = double>
class Vector;

// Forward declare the Matrix class for friend operators (ex. == and !=)
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E = double>
class Matrix;

// Declare operators for friendship
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
bool operator==(const Matrix<ROWS, COLS, E> & _m1, const Matrix<ROWS, COLS, E> & _m2);

template <std::size_t ROWS, std::size_t COLS, typename E>
bool operator!=(const Matrix<ROWS, COLS, E> & _m1, const Matrix<ROWS, COLS, E> & _m2);

template <std::size_t ROWS, std::size_t COLS, typename E>
std::ostream & operator<< (std::ostream & _stream, const Matrix<ROWS, COLS, E> & _matrix);


// Matrix class declaration
//--------------------------------------------------------------------------------------------
template
<std::size_t ROWS, std::size_t COLS, typename E>
class Matrix
{
    // Make every instantiation friend for inter-size opearations (transpose, minor)
    template <std::size_t, std::size_t, typename>
    friend class Matrix;
    
    // Allow vector class to access arr
    template <std::size_t, typename>
    friend class Vector;

private:
    E * arr;
    
    E & _get(std::size_t row, std::size_t column);
    const E & _get(std::size_t row, std::size_t column) const;

public:
    // Default, copy and move ctor, copy and move assingment operation, dtor
    //----------------------------------------------------------------------------------------
    Matrix();
    Matrix(const Matrix<ROWS, COLS, E> & _other);
    Matrix(Matrix<ROWS, COLS, E> && _other);
    Matrix & operator=(const Matrix<ROWS, COLS, E> & _other);
    Matrix & operator=(Matrix<ROWS, COLS, E> && _other);
    ~Matrix();
    
    // Other ctors
    //----------------------------------------------------------------------------------------
    Matrix(E _fill);
    Matrix(E **_arr);
    Matrix(E *_arr);
    Matrix(const std::array<E, ROWS * COLS> & _arr);

    // Friend * operator
    //----------------------------------------------------------------------------------------
    template <std::size_t ROWS1, std::size_t ROWS_COLS, std::size_t COLS2, typename _E>
    friend Matrix<ROWS1, COLS2, _E> operator*(const Matrix<ROWS1, ROWS_COLS, _E> & matrix1, const Matrix<ROWS_COLS, COLS2, _E> & matrix2);

    // Operators
    //----------------------------------------------------------------------------------------
    Matrix<ROWS, COLS, E> & operator+=(const Matrix<ROWS, COLS, E> & other);
    Matrix<ROWS, COLS, E> & operator/=(E x);
    Matrix<ROWS, COLS, E> & operator-=(const Matrix<ROWS, COLS, E> & other);
    Matrix<ROWS, COLS, E> & operator*=(const Matrix<COLS, COLS, E> & other);
    Matrix<ROWS, COLS, E> & operator*=(E x);
    Matrix<ROWS, COLS, E> operator-(void) const;

    E & operator()(std::size_t row, std::size_t column);
    const E & operator()(std::size_t row, std::size_t column) const;
    
    // Friendship declarations
    //----------------------------------------------------------------------------------------
    friend bool operator== <> (const Matrix & _m1, const Matrix & _m2);
    friend bool operator!= <> (const Matrix & _m1, const Matrix & _m2);
    friend std::ostream & operator<< <> (std::ostream & _stream, const Matrix & _matrix);

    // Other operations
    //----------------------------------------------------------------------------------------
    std::size_t getRows() const;
    std::size_t getCols() const;
    Matrix<COLS, ROWS, E> transpose() const;
    std::string toString() const;
    
    // Identity matrix
    //----------------------------------------------------------------------------------------
    template<std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    static typename std::enable_if<_ROWS == _COLS, Matrix<ROWS, COLS, E>>::type 
    identity();
    
    // Family of det methods
    //----------------------------------------------------------------------------------------
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && _ROWS == 1, E>::type
    det() const;
    
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && _ROWS == 2, E>::type
    det() const;
    
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && _ROWS == 3, E>::type
    det() const;
    
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && (_ROWS > 3), E>::type
    det() const;
    
    // Family of inverse methods
    //----------------------------------------------------------------------------------------
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && _ROWS == 1, Matrix<ROWS, COLS, E>>::type
    inverse() const;
    
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && _ROWS == 2, Matrix<ROWS, COLS, E>>::type
    inverse() const;
    
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && (_ROWS > 2), Matrix<ROWS, COLS, E>>::type
    inverse() const;
    
    // Family of matrix_minor methods
    //----------------------------------------------------------------------------------------
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && _ROWS == 2, E>::type
    matrix_minor(std::size_t _row, std::size_t _column) const;
    
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && (_ROWS > 2), E>::type
    matrix_minor(std::size_t _row, std::size_t _column) const;
    
    // Rotation matrices
    //----------------------------------------------------------------------------------------
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS, typename _E = E>
    static typename std::enable_if<_ROWS == 2 && _COLS == 2 && std::is_same<_E, double>::value, Matrix<ROWS, COLS, E>>::type
    rotation(double _a);
    
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS, typename _E = E>
    static typename std::enable_if<_ROWS == 3 && _COLS == 3 && std::is_same<_E, double>::value, Matrix<ROWS, COLS, E>>::type
    rotation(double _ax, double _ay, double _az);
};



#include "Matrix.tpp"


#endif  // _MATRIX_H
