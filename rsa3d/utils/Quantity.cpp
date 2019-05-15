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

// Returns Quantity computed from vector of samples
//----------------------------------------------------------------------------------------
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

// Friend stream insertion operator
//----------------------------------------------------------------------------------------
std::ostream &operator<<(std::ostream &stream, const Quantity &quantity)
{
    stream << quantity.value << " +- " << quantity.error;
    return stream;
}