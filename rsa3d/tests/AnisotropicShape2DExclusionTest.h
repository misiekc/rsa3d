//
// Created by PKua on 08.02.18.
//

#ifndef RSA3D_SPHEROCYLINDER2DINTTEST_H
#define RSA3D_SPHEROCYLINDER2DINTTEST_H

#include "../Vector.h"
#include "utility/ShapePairFactory.h"
#include "../AnisotropicShape2D.h"
#include <memory>

namespace as2d_extest
{
    struct Point {
        Vector<2> position;
        std::string color;
    };

    struct Context {
        enum Mode {
            OVERLAP,
            POINT_INSIDE
        };

        std::string particle;
        std::string attributes;
        Mode mode{};
        Vector<2> firstPos;
        double firstAngle{};
        double box_width{};
        double box_height{};
        double points{};
        bool append{};

        // For Mode::OVERLAP
        double secondAngle{};

        // For Mode::POINT_INSIDE
        double fromAngle{};
        double toAngle{};

        void load(const std::string & _filename);
        void printInfo() const;
        Point randomPoint(RND * rnd) const;
    };

    void perform_inttest(const Context &_context, std::vector<Point> &_result);
    void perform_pitest(const Context &_context, std::vector<Point> &_result);
    void results_to_wolfram(const Context &_context, const std::vector<as2d_extest::Point> &_results,
                            std::ostream &_out);
    int main(int argc, char **argv);

    Point randomPoint(const Context &_context, RND &rnd);
    Point findMaxXRedPoint(const std::vector<Point> &_results);

    std::unique_ptr<AnisotropicShape2D> generateShape(double angle, const Vector<2> &position);
}

#endif //RSA3D_SPHEROCYLINDER2DINTTEST_H
