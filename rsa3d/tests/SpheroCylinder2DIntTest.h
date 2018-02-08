//
// Created by PKua on 08.02.18.
//

#ifndef RSA3D_SPHEROCYLINDER2DINTTEST_H
#define RSA3D_SPHEROCYLINDER2DINTTEST_H

#include "../Vector.h"
#include "utility/ShapePairFactory.h"

namespace spheroc_inttest
{
    struct Point {
        Vector<2> position;
        std::string color;
    };

    struct Context {
        Vector<2> firstPos;
        double firstAngle{};
        double secondAngle{};
        double ratio{};
        double box_width{};
        double box_height{};
        double points{};

        void load(const std::string & _filename);
    };

    void perform(const Context &_context, std::vector<Point> &_result);
    void print_results(const Context &_context, const std::vector<spheroc_inttest::Point> &_results,
                           std::ostream &_out);
    int main(int argc, char **argv);
}

#endif //RSA3D_SPHEROCYLINDER2DINTTEST_H
