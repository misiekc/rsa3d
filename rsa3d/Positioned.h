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

template <unsigned short SPATIAL_DIMENSION>
class Positioned {
private:
    std::array<double, SPATIAL_DIMENSION> position;

protected:

	/**
     * Sets new position. Derived shape classes have to override this method when they want to keep track of
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
	 * Returns position of a Positioned.
	 * @return position of a Positioned
	 */
	double* getPosition() const;

    /**
     * Translates positioned by a given vector @a v.
     * @param v
     */
	void translate(double *v);
};

#include "Positioned.tpp"

#endif /* POSITIONED_H_ */
