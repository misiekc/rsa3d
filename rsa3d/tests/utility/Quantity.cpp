//
// Created by PKua on 25.12.17.
//

#include "../CuboidSpeedTest.h"
#include <fstream>
#include <bits/unique_ptr.h>
#include <iomanip>
#include <chrono>

// Returns Quantity computed from vector of samples
//----------------------------------------------------------------------------------------
Quantity Quantity::fromSamples(const std::vector<double> & _samples)
{
    if (_samples.empty())
        return Quantity();
    else if (_samples.size() == 1)
        return Quantity(_samples[0], 0);

    Quantity result;
    double sum, dev_sum;

    sum = accumulate(_samples.begin(), _samples.end(), 0);
    result.value = sum / _samples.size();
    dev_sum = accumulate(_samples.begin(), _samples.end(), 0,
        [&](double s, double next) {
            s += pow(result.value - next, 2);
            return s;
        });
   result.error = dev_sum / _samples.size() / (_samples.size() - 1);
   result.error = sqrt(result.error);
   return result;
}

// Friend stream insertion operator
//----------------------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & _stream, const Quantity & _quantity)
{
    _stream << _quantity.value << " +- " << _quantity.error;
    return _stream;
}