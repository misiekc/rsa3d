//-----------------------------------------------------------------------------------------------------------------------------
// Class representing matrices and operations on them. Template definitions
//-----------------------------------------------------------------------------------------------------------------------------
// (C)PKua 2017
//-----------------------------------------------------------------------------------------------------------------------------


// Default ctor generating zero-filled matrix
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E>::Matrix() : Matrix<ROWS, COLS, E>(E(0))
{

}


// Copy ctor
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E>::Matrix(const Matrix<ROWS, COLS, E> & _other)
{
    // Alloc an array and copy data
    arr = new E[ROWS * COLS];
    std::copy(_other.arr, _other.arr + (ROWS * COLS), arr);
}


// Move ctor
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E>::Matrix(Matrix<ROWS, COLS, E> && _other) :
    arr(_other.arr)
{
    _other.arr = nullptr;
}


// Copy assingment operator. Returns this reference for operations such as
// (matrix = getmatrix()).getWidth()
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E> & Matrix<ROWS, COLS, E>::operator=(const Matrix<ROWS, COLS, E> & _other)
{
    // Self assingment, skip
    if (this == &_other)
        return *this;

    // Copy data
    std::copy(_other.arr, _other.arr + (ROWS * COLS), arr);

    return *this;
}


// Move assignment operator. Returns this reference for operations such as
// (matrix = getmatrix()).getWidth()
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E> & Matrix<ROWS, COLS, E>::operator=(Matrix<ROWS, COLS, E> && _other)
{
    delete[] arr;
    arr = _other.arr;
    _other.arr = nullptr;

    return *this;
}


// Dtor
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E>::~Matrix()
{
    delete[] arr;
}


// Ctor creating matrix of given size. Fills it with given value
//--------------------------------------------------------------------------------------------
// _fill - value to be filled with
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E>::Matrix(E _fill)
{
    // Alloc an array and fill
    arr = new E[ROWS * COLS];
    std::fill(arr, arr + (ROWS * COLS), _fill);
}


// Ctor creating matrix from given two-dimentional array, where the element from the i-th row
// and the j-th column is in arr[i-1][j-1]. Data from array are copied
//--------------------------------------------------------------------------------------------
// _arr - array with matrix elements
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E>::Matrix(E **_arr)
{
    // Alloc an array and copy elements from given
    arr = new E[ROWS * COLS];
    std::size_t arr_index = 0;
    for (std::size_t i = 0; i < ROWS; i++)
        for (std::size_t j = 0; j < COLS; j++)
            arr[arr_index++] = _arr[i][j];
}


// Ctor creating matrix from given one-dimentional array. When desired matrix is
// / 1 2 3 |
// | 4 5 6 |
// \ 7 8 9 |
// then the array should be: {1, 2, 3, 4, 5, 6, 7, 8, 9}
//--------------------------------------------------------------------------------------------
// _arr - array with matrix elements
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E>::Matrix(E *_arr)
{
    // Alloc an array and copy elements from given
    arr = new E[ROWS * COLS];
    std::copy(_arr, _arr + (ROWS * COLS), arr);
}


// Ctor creating matrix from std::array. When desired matrix is
// / 1 2 3 |
// | 4 5 6 |
// \ 7 8 9 |
// then the array should be: {1, 2, 3, 4, 5, 6, 7, 8, 9}
//--------------------------------------------------------------------------------------------
// _arr - initializer lsit with matrix elements
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E>::Matrix(const std::array<E, ROWS * COLS> & _arr)
{
    // Alloc an array and copy elements from given
    arr = new E[ROWS * COLS];
    std::copy(_arr.begin(), _arr.begin() + (ROWS * COLS), arr);
}


// Addition of 2 matrices
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
inline Matrix<ROWS, COLS, E> operator+(Matrix<ROWS, COLS, E> matrix1, const Matrix<ROWS, COLS, E> & matrix2)
{
    return (matrix1 += matrix2);
}


// Subtraction of 2 matrices
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
inline Matrix<ROWS, COLS, E> operator-(Matrix<ROWS, COLS, E> matrix1, const Matrix<ROWS, COLS, E> & matrix2)
{
    return (matrix1 -= matrix2);
}


// Multiplying matrix by scalar
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
inline Matrix<ROWS, COLS, E> operator*(Matrix<ROWS, COLS, E> matrix, E x)
{
    return (matrix *= x);
}


template <std::size_t ROWS, std::size_t COLS, typename E>
inline Matrix<ROWS, COLS, E> operator*(E x, Matrix<ROWS, COLS, E> matrix)
{
    return (matrix *= x);
}


// Dividing matrix by scalar
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
inline Matrix<ROWS, COLS, E> operator/(Matrix<ROWS, COLS, E> matrix, E x)
{
    return (matrix /= x);
}


// Multiplication of 2 matrices
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS1, std::size_t ROWS_COLS, std::size_t COLS2, typename E>
Matrix<ROWS1, COLS2, E> operator*(const Matrix<ROWS1, ROWS_COLS, E> & matrix1, const Matrix<ROWS_COLS, COLS2, E> & matrix2)
{
    Matrix<ROWS1, COLS2, E> ret;
    for (std::size_t i = 0; i < ROWS1; i++)
        for (std::size_t j = 0; j < COLS2; j++)
            for (std::size_t k = 0; k < ROWS_COLS; k++)
                ret._get(i, j) += matrix1._get(i, k) * matrix2._get(k, j);
    return ret;
}


// Equality operator of 2 matrices
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
bool operator==(const Matrix<ROWS, COLS, E> & matrix1, const Matrix<ROWS, COLS, E> & matrix2)
{
    for (std::size_t i = 0; i < ROWS * COLS; i++)
        if (matrix1.arr[i] != matrix2.arr[i])
            return false;
    return true;
}


// Inequality operator of 2 matrices
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
inline bool operator!=(const Matrix<ROWS, COLS, E> & matrix1, const Matrix<ROWS, COLS, E> & matrix2)
{
    return !(matrix1 == matrix2);
}


// Stream insertion operator
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
std::ostream & operator<< (std::ostream & _stream, const Matrix<ROWS, COLS, E> & _matrix)
{
    for (std::size_t i = 0; i < ROWS; i++) {
        for (std::size_t j = 0; j < COLS; j++) {
            _stream << _matrix._get(i, j) << " ";
        }
        _stream << std::endl;
    }

    return _stream;
}


// Addition assignment operator
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E> & Matrix<ROWS, COLS, E>::operator+=(const Matrix<ROWS, COLS, E> & other)
{
    for (std::size_t i = 0; i < ROWS * COLS; i++)
        arr[i] += other.arr[i];
    return *this;
}


// Matrix multiplication assignment operator
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
inline Matrix<ROWS, COLS, E> & Matrix<ROWS, COLS, E>::operator*=(const Matrix<COLS, COLS, E> & other)
{
    *this = std::move(*this * other);
    return *this;
}


// Matrix multiplication by scalar assignment operator
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E> & Matrix<ROWS, COLS, E>::operator*=(E x)
{
    if (x == 0) {
        std::fill(arr, arr + (ROWS * COLS), 0);
    } else if (x != 1) {
        for (std::size_t i = 0; i < ROWS * COLS; i++)
            arr[i] *= x;
    }
    return *this;
}


// Matrix division by scalar assignment operator
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E> & Matrix<ROWS, COLS, E>::operator/=(E x)
{
    if (x != 1)
        for (std::size_t i = 0; i < ROWS * COLS; i++)
            arr[i] /= x;
    return *this;
}


// Subtraction assignment operator
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E> & Matrix<ROWS, COLS, E>::operator-=(const Matrix & other)
{
    for (std::size_t i = 0; i < ROWS * COLS; i++)
        arr[i] -= other.arr[i];
    return *this;
}


// Unary minus operator
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E> Matrix<ROWS, COLS, E>::operator-(void) const
{
    Matrix<ROWS, COLS, E> ret;
    for (std::size_t i = 0; i < ROWS * COLS; i++)
        ret.arr[i] = -arr[i];
    return ret;
}

// R/W access oparator. In order to get the element from i-th row and j-th column, one should
// use matrix(i-1, j-1). Check if indexes are in bounds and throws std::invalid_argument
// exception if needed
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
inline E & Matrix<ROWS, COLS, E>::operator()(std::size_t _row, std::size_t _column)
{
    if (_row >= ROWS)       throw std::invalid_argument("row index out of bounds");
    if (_column >= COLS)    throw std::invalid_argument("column index out of bounds");

    return _get(_row, _column);
}

// Read-only access oparator. In order to get the element from i-th row and j-th column, one
// should use matrix(i-1, j-1). Check if indexes are in bounds and throws
// std::invalid_argument exception if needed
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
inline const E & Matrix<ROWS, COLS, E>::operator()(std::size_t _row, std::size_t _column) const
{
    if (_row >= ROWS)       throw std::invalid_argument("row index out of bounds");
    if (_column >= COLS)    throw std::invalid_argument("column index out of bounds");

    return _get(_row, _column);
}


// Returns transposition of the matrix
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<COLS, ROWS, E> Matrix<ROWS, COLS, E>::transpose() const
{
    Matrix<COLS, ROWS, E> matrix;
    for (std::size_t i = 0; i < ROWS; i++)
        for (std::size_t j = 0; j < COLS; j++)
            matrix._get(j, i) = _get(i, j);

    return matrix;
}


// Returns inversion of the matrix. For 1x1 matrix
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
template <std::size_t _ROWS, std::size_t _COLS>
typename std::enable_if<_ROWS == _COLS && _ROWS == 1, Matrix<ROWS, COLS, E>>::type
Matrix<ROWS, COLS, E>::inverse() const
{
    if (arr[0] == E(0))  throw std::runtime_error("cannot invert, det == 0");
    return Matrix<1, 1, E> (1 / arr[0]);
}


// Returns inversion of the matrix. For 2x2 matrix
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
template <std::size_t _ROWS, std::size_t _COLS>
typename std::enable_if<_ROWS == _COLS && _ROWS == 2, Matrix<ROWS, COLS, E>>::type
Matrix<ROWS, COLS, E>::inverse() const
{
    E determinant = det();
    Matrix<2, 2, E> out;
    if (determinant == E(0))  throw std::runtime_error("cannot invert, det == 0");
    
    out.arr[0] = arr[3] / determinant;
    out.arr[1] = -arr[1] / determinant;
    out.arr[2] = -arr[2] / determinant;
    out.arr[3] = arr[0] / determinant;
    return out;
}


// Returns inversion of the matrix. For 3x3 and bigger matrices
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
template <std::size_t _ROWS, std::size_t _COLS>
typename std::enable_if<_ROWS == _COLS && (_ROWS > 2), Matrix<ROWS, COLS, E>>::type
Matrix<ROWS, COLS, E>::inverse() const
{
    E determinant = det();
    Matrix<ROWS, COLS, E> out;
    if (determinant == E(0))  throw std::runtime_error("cannot invert, det == 0");

    for (std::size_t i = 0; i < ROWS; i++)
        for (std::size_t j = 0; j < COLS; j++)
            out._get(i, j) = ((i + j) % 2 == 0) ? (matrix_minor(j, i) / determinant) : (-matrix_minor(j, i) / determinant);
    return out;
}


// Returns identity matrix of template parameters size
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
template <std::size_t _ROWS, std::size_t _COLS>
typename std::enable_if<_ROWS == _COLS, Matrix<ROWS, COLS, E>>::type 
Matrix<ROWS, COLS, E>::identity()
{
    Matrix<ROWS, COLS, E> matrix(E(0));
    for (std::size_t i = 0; i < COLS; i++)
        matrix._get(i, i) = E(1);
    return matrix;
}


// Calculates minor of matrix created by removing (_row + 1) row and (_column + 1) column.
// For 2x2 matrix
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
template <std::size_t _ROWS, std::size_t _COLS>
typename std::enable_if<_ROWS == _COLS && _ROWS == 2, E>::type
Matrix<ROWS, COLS, E>::matrix_minor(std::size_t _row, std::size_t _column) const
{   
    return arr[3 ^ (_row << 1 | _column)];
}


// Calculates minor of matrix created by removing (_row + 1) row and (_column + 1) column.
// For 3x3 and bigger matrices
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
template <std::size_t _ROWS, std::size_t _COLS>
typename std::enable_if<_ROWS == _COLS && (_ROWS > 2), E>::type
Matrix<ROWS, COLS, E>::matrix_minor(std::size_t _row, std::size_t _column) const
{   
    Matrix<ROWS - 1, COLS - 1, E> minor_mx;

    if (_row > 0) {
        // Copy left-top part (if exists)
        if (_column > 0)
            for (std::size_t i = 0; i < _row; i++)
                for (std::size_t j = 0; j < _column; j++)
                    minor_mx._get(i, j) = _get(i, j);

        // Copy right-top part (if exists)
        if (_column < COLS - 1)
            for (std::size_t i = 0; i < _row; i++)
                for (std::size_t j = _column + 1; j < COLS; j++)
                    minor_mx._get(i, j - 1) = _get(i, j);
    }

    if (_row < ROWS - 1) {
        // Copy left-bottom part (if exists)
        if (_column > 0)
            for (std::size_t i = _row + 1; i < ROWS; i++)
                for (std::size_t j = 0; j < _column; j++)
                    minor_mx._get(i - 1, j) = _get(i, j);

        // Copy right-bottom part (if exists)
        if (_column < COLS - 1)
            for (std::size_t i = _row + 1; i < ROWS; i++)
                for (std::size_t j = _column + 1; j < COLS; j++)
                    minor_mx._get(i - 1, j - 1) = _get(i, j);
    }

    return minor_mx.det();
}


// Returns determinant of the matrix. For 1x1 matrix
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
template <std::size_t _ROWS, std::size_t _COLS>
typename std::enable_if<_ROWS == _COLS && _ROWS == 1, E>::type
Matrix<ROWS, COLS, E>::det() const
{
    return arr[0];
}


// Returns determinant of the matrix. For 2x2 matrix
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
template <std::size_t _ROWS, std::size_t _COLS>
typename std::enable_if<_ROWS == _COLS && _ROWS == 2, E>::type
Matrix<ROWS, COLS, E>::det() const
{
    return arr[0] * arr[3] - arr[1] * arr[2];
}


// Returns determinant of the matrix. For 3x3 matrix
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
template <std::size_t _ROWS, std::size_t _COLS>
typename std::enable_if<_ROWS == _COLS && _ROWS == 3, E>::type
Matrix<ROWS, COLS, E>::det() const
{
    return arr[0] * arr[4] * arr[8] + arr[1] * arr[5] * arr[6] + arr[2] * arr[3] * arr[7]
        - arr[2] * arr[4] * arr[6] - arr[1] * arr[3] * arr[8] - arr[0] * arr[5] * arr[7];
}


// Returns determinant of the matrix. For 4x4 and bigger matrices
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
template <std::size_t _ROWS, std::size_t _COLS>
typename std::enable_if<_ROWS == _COLS && (_ROWS > 3), E>::type
Matrix<ROWS, COLS, E>::det() const
{
    // Laplace expansion
    E determinant = E(0.0);
    bool even_perm = true;
    for (std::size_t i = 0; i < ROWS; i++) {
        if (even_perm)  determinant += _get(0, i) * matrix_minor(0, i);
        else            determinant -= _get(0, i) * matrix_minor(0, i);
        even_perm = !even_perm;
    }
    return determinant;
}


// Generates two-dimentional rotate matrix
//--------------------------------------------------------------------------------------------
// _a - counter-clockwise rotate angle
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
template <std::size_t _ROWS, std::size_t _COLS, typename _E>
typename std::enable_if<_ROWS == 2 && _COLS == 2 && std::is_same<_E, double>::value, Matrix<ROWS, COLS, E>>::type
Matrix<ROWS, COLS, E>::rotation(double _a)
{
    double sin_a = sin(_a);
    double cos_a = cos(_a);
    Matrix<2, 2, double> matrix;
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
template <std::size_t ROWS, std::size_t COLS, typename E>
template <std::size_t _ROWS, std::size_t _COLS, typename _E>
typename std::enable_if<_ROWS == 3 && _COLS == 3 && std::is_same<_E, double>::value, Matrix<ROWS, COLS, E>>::type
Matrix<ROWS, COLS, E>::rotation(double _ax, double _ay, double _az)
{
    double sin_ax = sin(_ax);
    double sin_ay = sin(_ay);
    double sin_az = sin(_az);
    double cos_ax = cos(_ax);
    double cos_ay = cos(_ay);
    double cos_az = cos(_az);

    Matrix<3, 3, double> matrix;
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


// Depicts the matrix in the form of a string
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
std::string Matrix<ROWS, COLS, E>::toString() const
{
    std::stringstream sstream;
    sstream << *this;
    return sstream.str();
}

// Inline function returning number of rows in matrix
//-------------------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
inline std::size_t Matrix<ROWS, COLS, E>::getRows() const
{
    return ROWS;
}

// Inline function returning number of columns in matrix
//-------------------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
inline std::size_t Matrix<ROWS, COLS, E>::getCols() const
{
    return COLS;
}

// Private inline function returning the element from the array of given address. Doesn't perform
// arguments validation
//-------------------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
inline E & Matrix<ROWS, COLS, E>::_get(std::size_t _row, std::size_t _col)
{
    return arr[_row * COLS + _col];
}

// Private inline function returning immutable element from the array of given address. Doesn't perform
// arguments validation
//-------------------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E>
inline const E & Matrix<ROWS, COLS, E>::_get(std::size_t _row, std::size_t _col) const
{
    return arr[_row * COLS + _col];
}