//
// Created by PKua on 25.12.17.
//

#include "Quantity.h"
#include <fstream>
#include <bits/unique_ptr.h>
#include <iomanip>
#include <chrono>
#include <vector>
#include <numeric>
#include <cmath>


Quantity Quantity::fromSamples(const std::vector<double> &samples)
{
    if (samples.empty())
        return Quantity();
    else if (samples.size() == 1)
        return Quantity(samples[0], 0);

    Quantity result;

    result.value = std::accumulate(samples.begin(), samples.end(), 0.0) / samples.size();

    auto squareDev = [&](double s, double next) { return s + std::pow(result.value - next, 2); } ;
    result.error = std::accumulate(samples.begin(), samples.end(), 0.0, squareDev);
    result.error = result.error / samples.size() / (samples.size() - 1);
    result.error = std::sqrt(result.error);

    return result;
}

std::ostream &operator<<(std::ostream &stream, const Quantity &quantity)
{
    std::string value, error;
    if (!quantity.significantDigitsBasedOnError) {
        value = std::to_string(quantity.value);
        error = std::to_string(quantity.error);
    } else {

        /* Here goes the implementation */

    }

    stream << value << quantity.getSeparatorString() << error;
    return stream;
}

std::string Quantity::getSeparatorString() const {
    switch (this->separator) {
        case TABULATOR:
            return "\t";
        case PLUS_MINUS:
            return " +- ";
    }
    return "";
}