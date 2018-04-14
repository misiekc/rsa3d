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

/**
 * @brief An object located in @a SPATIAL_DIMENSION dimensional space.
 *
 * The class abstracts any object which position can be described by @a SPATIAL_DIMENSION numbers. Its position can
 * be accessed or translated by a given vector.
 * @tparam SPATIAL_DIMENSION number of dimensions of a space in which positioned lives
 */
template <unsigned short SPATIAL_DIMENSION>
class Positioned {
private:
    std::array<double, SPATIAL_DIMENSION> position;

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
    virtual void setPosition(const double *position);

public:
	virtual ~Positioned() = default;

	/**
	 * @brief Returns position of a Positioned.
	 * @return position of a Positioned
	 */
	double* getPosition() const;

    /**
     * @brief Translates positioned by a given vector @a v.
     * @param v a vector to translate by
     */
	void translate(double *v);
};

#include "Positioned.tpp"

#endif /* POSITIONED_H_ */
