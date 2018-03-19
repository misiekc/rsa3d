//
// Created by PKua on 25.12.17.
//

#ifndef RSA3D_QUANTITY_H
#define RSA3D_QUANTITY_H

#include <vector>
#include <iosfwd>

// Struct representing quantity with uncertainity
struct Quantity {
    double value = 0;
    double error = 0;

    Quantity() = default;

    Quantity(double _value, double _error) : value(_value), error(_error)
    { }

    static Quantity fromSamples(const std::vector<double> & _samples);
    friend std::ostream & operator<<(std::ostream & _stream, const Quantity & _quantity);
};

#endif //RSA3D_QUANTITY_H
