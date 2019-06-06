//
// Created by PKua on 25.12.17.
//

#ifndef RSA3D_QUANTITY_H
#define RSA3D_QUANTITY_H

#include <vector>
#include <iosfwd>

/**
 * @brief A simple class representing physical quantity with error.
 *
 * The class also includes simple print routines, where the number of significant digits can be based on error.
 */
class Quantity {

public:
    /**
     * @brief Enumeration representing what sould be used to separate value from error when printing.
     */
    enum Separator {
        TABULATOR,
        PLUS_MINUS
    };

    /**
     * @brief The value of a quantity.
     */
    double value = 0;

    /**
     * @brief The error of a quantiy.
     */
    double error = 0;

    /**
     * @brief What should be used to separate value from error when printing.
     */
    Separator separator = TABULATOR;

    /**
     * @brief It set true, the number of significant digits will be based on error. If false, default double printing.
     */
    bool significantDigitsBasedOnError = false;

    Quantity() = default;
    Quantity(double value, double error) : value(value), error(error) { }

    /**
     * @brief Creates the quantity as a mean of samples, with estimated error of the mean.
     * @param samples A vector of doubles to calculate the quantity from
     * @return the quantity from samples
     */
    static Quantity fromSamples(const std::vector<double> &samples);

    /**
     * @brief Prints the quantity on given @a std::ostream The behaviour can be manipulated via Quantity::separator and
     * Quantity::significantDigitsBasedOnError fields.
     * @param stream the @a std::ostream to print onto
     * @param quantity the quantity printed
     * @return the reference to `stream`
     */
    friend std::ostream &operator<<(std::ostream &stream, const Quantity &quantity);

private:
    std::string getSeparatorString() const;
};

#endif //RSA3D_QUANTITY_H
