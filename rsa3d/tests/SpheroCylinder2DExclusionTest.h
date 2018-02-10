//
// Created by PKua on 08.02.18.
//

#ifndef RSA3D_SPHEROCYLINDER2DINTTEST_H
#define RSA3D_SPHEROCYLINDER2DINTTEST_H

#include "../Vector.h"
#include "utility/ShapePairFactory.h"

namespace spheroc_extest
{
    struct Point {
        Vector<2> position;
        std::string color;
    };

    struct Context {
        Vector<2> firstPos;
        double firstAngle{};
        double fromAngle{};
        double toAngle{};
        double ratio{};
        double box_width{};
        double box_height{};
        double points{};
        bool append{};

        void load(const std::string & _filename);
    };

    void perform_inttest(const Context &_context, std::vector<Point> &_result);
    void perform_pitest(const Context &_context, std::vector<Point> &_result);
    void results_to_wolfram(const Context &_context, const std::vector<spheroc_extest::Point> &_results,
                            std::ostream &_out);
    int main(int argc, char **argv);
}

#endif //RSA3D_SPHEROCYLINDER2DINTTEST_H
