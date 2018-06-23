/*
 * Positioned.h
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#ifndef POSITIONED_H_
#define POSITIONED_H_

// #include "Vector.h"

#include <array>
#include "Vector.h"

/**
 * @brief An object located in @a SPATIAL_DIMENSION dimensional space.
 *
 * The class abstracts any object which position can be described by @a SPATIAL_DIMENSION numbers. Its position can
 * be accessed or translated by a given vector.
 * @tparam SPATIAL_DIMENSION number of dimensions of a space in which positioned lives
 */
template <unsigned short SPATIAL_DIMENSION>
class Positioned {
	static_assert(SPATIAL_DIMENSION > 0, "SPATIAL_DIMENTION == 0");

private:
    Vector<SPATIAL_DIMENSION> position;

protected:

	/**
     * @brief Sets new position.
     *
     * Derived shape classes have to override this method when they want to keep track of
     * Positioned position. One would then typically write:
     * \code
     * void Derived::setPosition(const double *position) {
     *     Positioned::setPosition(position);
     *     // invoke getPosition and for example compute vertices
     * }
     * \endcode
     * @param position new position
     */
    virtual void setPosition(const Vector<SPATIAL_DIMENSION> &position);

public:

    static int offset[(1 << SPATIAL_DIMENSION)][SPATIAL_DIMENSION];

	virtual ~Positioned() = default;

	/**
	 * @brief Returns position of a Positioned.
	 * @return position of a Positioned
	 */
	const Vector<SPATIAL_DIMENSION> &getPosition() const;

    /**
     * @brief Translates positioned by a given vector @a v.
     * @param v a vector to translate by
     */
	void translate(const Vector<SPATIAL_DIMENSION> &position);
};

using RSAPositioned = Positioned<RSA_SPATIAL_DIMENSION>;

#include "Positioned.tpp"

#endif /* POSITIONED_H_ */
