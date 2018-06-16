//
// PKua on 08.02.18.
//

#include <fstream>
#include <bits/unique_ptr.h>
#include <memory>
#include "AnisotropicShape2DExclusionTest.h"
#include "../rsa3d/Config.h"
#include "../rsa3d/ShapeFactory.h"
#include "../rsa3d/Utils.h"
#include "utility/MockBC.h"
#include "utility/ShapeGenerators.h"

namespace
{
    const char *const FIRST_SHAPE_COLOR = "Black";
    const char *const SECOND_SHAPE_COLOR = "Blue";
    const char *const INSIDE_COLOR = "Orange";
    const char *const OUTSIDE_COLOR = "Green";
    const char *const CIRCLE_FIX_COLOR = "Purple";
    const char *const FIRST_ZONE_COLOR = "Red";
    const char *const SECOND_ZONE_COLOR = "Yellow";

    struct ColoredPoint {
        Vector<2> position;
        std::string color;
    };

    using Result = std::vector<ColoredPoint>;

    /* Settings for test loaded from file and some utilities*/
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
        bool wolframSupported{};

        // For Mode::OVERLAP
        double secondAngle{};

        // For Mode::POINT_INSIDE
        double fromAngle{};
        double toAngle{};

        void load(const std::string &filename);
        void printInfo(std::ostream &out) const;
        ColoredPoint randomPoint(RND *rnd) const;
    };

    /* Loads context from a given file */
    void Context::load(const std::string &filename) {
        std::ifstream input(filename);
        if (!input)     die("Error opening " + std::string(filename) + " file to read");

        auto config = std::unique_ptr<Config>(Config::parse(input));
        this->particle = config->getString("particle");
        this->attributes = config->getString("attr");
        ShapeFactory::initShapeClass(this->particle, this->attributes);

        this->firstPos = Vector<2>{{config->getDouble("first_x"), config->getDouble("first_y")}};
        this->firstAngle = config->getDouble("first_angle") * M_PI / 180;

        std::string modeStr = config->getString("mode");
        if (modeStr == "overlap") {
            this->mode = Mode::OVERLAP;
            this->secondAngle = config->getDouble("second_angle") * M_PI / 180;
        } else if (modeStr == "point_inside") {
            this->mode = Mode::POINT_INSIDE;
            this->fromAngle = config->getDouble("from_angle") * M_PI / 180;
            this->toAngle = config->getDouble("to_angle") * M_PI / 180;
        } else {
            die("Unknown test mode: " + modeStr);
        }

        this->box_width = config->getDouble("box_width");
        this->box_height = config->getDouble("box_height");
        this->points = config->getDouble("points");
        this->append = config->getString("append") == "true";
        this->wolframSupported = !generate_shape<2, 1>(Vector<2>(), {{0}})->toWolfram().empty();
        if (!this->wolframSupported)
            std::cout << "[WARNING] " << this->particle << " has no toWolfram(). Shapes will not be drawn." << std::endl;
        std::cout << "[INFO] Context loaded from " << filename << std::endl;
    }

    /* Performs inttest based on context and returns result */
    Result perform_inttest(const Context &context, std::ostream &out) {
        RND rnd;
        MockBC mockBC;
        Result result;

        context.printInfo(out);
        auto first = generate_shape<2, 1>(context.firstPos, {{context.firstAngle}});
        for (std::size_t i = 0; i < context.points; i++) {
            ColoredPoint point = context.randomPoint(&rnd);
            auto second = generate_shape<2, 1>(point.position, {{context.secondAngle}});

            if (first->overlap(&mockBC, second.get()))
                point.color = INSIDE_COLOR;
            else
                point.color = OUTSIDE_COLOR;
            result.push_back(point);
        }
        return result;
    }

    /* Uses pointInside method when dealing with ConvexShape or voxelInside with zero spatial size for normal Shape */
    bool point_inside(const Shape<2, 1> &shape, const Vector<2> &point, double angleFrom, double angleTo) {
        MockBC bc;
        double pointArr[2];
        point.copyToArray(pointArr);

        try {
            auto &convexShape = dynamic_cast<const ConvexShape<2, 1>&>(shape);
            return convexShape.pointInside(&bc, pointArr, {{angleFrom}}, angleTo - angleFrom);
        } catch (std::bad_cast&) {
            return shape.voxelInside(&bc, pointArr, {{angleFrom}}, 0, angleTo - angleFrom);
        }
    }

    /* Performs pitest based on context and returns result */
    Result perform_pitest(const Context &context, std::ostream &out) {
        MockBC mockBC;
        RND rnd;
        Result result;

        context.printInfo(out);
        auto shape = generate_shape<2, 1>(context.firstPos, {{context.firstAngle}});
        for (std::size_t i = 0; i < context.points; i++) {
            ColoredPoint point = context.randomPoint(&rnd);
            if (point_inside(*shape, point.position, context.fromAngle, context.toAngle)) {
                point.color = INSIDE_COLOR;
            } else {
                auto angleFromShape = generate_shape<2, 1>(point.position, {{context.fromAngle}});
                auto angleToShape = generate_shape<2, 1>(point.position, {{context.toAngle}});
                auto withinAngleFromZone = shape->overlap(&mockBC, angleFromShape.get());
                auto withinAngleToZone = shape->overlap(&mockBC, angleToShape.get());

                if (withinAngleFromZone && withinAngleToZone)
                    point.color = CIRCLE_FIX_COLOR;
                else if (withinAngleFromZone)
                    point.color = FIRST_ZONE_COLOR;
                else if (withinAngleToZone)
                    point.color = SECOND_ZONE_COLOR;
                else
                    point.color = OUTSIDE_COLOR;
            }
            result.push_back(point);
        }

        return result;
    }

    /* Checks in context what test to perform and executes it*/
    std::vector<ColoredPoint> perform_test(const Context &context, std::ostream &out) {
        if (context.mode == Context::Mode::OVERLAP)
            return perform_inttest(context, out);
        else  // Context::Mode::POINT_INSIDE
            return perform_pitest(context, out);
    }

    ColoredPoint Context::randomPoint(RND * rnd) const {
        ColoredPoint point;
        Vector<2> secondTranslation{{
            (rnd->nextValue() - 0.5) * box_width,
            (rnd->nextValue() - 0.5) * box_height
        }};
        point.position = firstPos + secondTranslation;
        return point;
    }

    /* Prints info about current context */
    void Context::printInfo(std::ostream &out) const {
        if (this->mode == Mode::OVERLAP)
            out << "[OVERLAP TEST] ";
        else    // Mode::POINT_INSIDE
            out << "[POINT INSIDE TEST] ";

        auto first = generate_shape<2, 1>(this->firstPos, {{this->firstAngle}});
        out << "Generating " << this->points << " points from " << this->box_width << " x ";
        out << this->box_height << " box for:" << std::endl;
        out << first->toString() << std::endl;

        if (this->mode == Mode::OVERLAP)
            out << "Second angle: " << this->secondAngle * 180 / M_PI << std::endl;
        else    // Mode::POINT_INSIDE
            out << "Angle range: [" << this->fromAngle * 180 / M_PI << ", " << this->toAngle * 180 / M_PI << "]" << std::endl;
    }

    /* Finds a red point in the result which X coordinate is the greatest - used for drawing the second shape*/
    ColoredPoint find_max_x_red_point(const Result &result) {
        auto comparator = [](const ColoredPoint & p1, const ColoredPoint & p2) {
            if (p1.color != INSIDE_COLOR && p2.color == INSIDE_COLOR)
                return true;
            else if (p1.color == INSIDE_COLOR && p2.color != INSIDE_COLOR)
                return false;
            else
                return p1.position[0] < p2.position[0];
        };

        return *max_element(result.begin(), result.end(), comparator);
    }

    /* Converts results to executable woflram notebook drawing all stuff */
    void results_to_wolfram(const Context &context, const Result &result, std::ostream &out) {
        if (context.append)
            out << std::endl << std::endl;

        out << std::fixed;
        out << "points = {" << std::endl;
        for (ColoredPoint point : Result(result.begin(), result.end() - 1))
            out << "    {" << point.color << ", Point[" << point.position << "]}," << std::endl;
        out << "    {" << result.back().color << ", Point[" << result.back().position << "]}};" << std::endl;

        auto first = generate_shape<2, 1>(context.firstPos, {{context.firstAngle}});
        if (!context.wolframSupported) {
            out << std::endl << "Graphics[points]";
        } else {
            out << "shape1 = " << first->toWolfram() << ";" << std::endl;
            if (context.mode == Context::Mode::OVERLAP) {
                auto second = generate_shape<2, 1>(find_max_x_red_point(result).position, {{context.secondAngle}});
                out << "shape2 = " << second->toWolfram() << ";" << std::endl << std::endl;
                out << "Graphics[Join[{{" << SECOND_SHAPE_COLOR << ", shape2}, {" << FIRST_SHAPE_COLOR << ", shape1}}, points]]";
            } else {  // Context::Mode::POINT_INSIDE
                out << std::endl << "Graphics[Join[{{" << FIRST_SHAPE_COLOR << ", shape1}}, points]]";
            }
        }
    }
}

int as2d_extest::main(int argc, char **argv) {
    if (argc < 3)
        die("Usage: ./rsa_test as2d_extest <input> <output = out.nb>");

    Context context;
    context.load(argv[2]);

    Result result = perform_test(context, std::cout);
    std::string filename = (argc == 4 ? argv[3] : "out.nb");
    std::ofstream output(filename, context.append ? std::ofstream::app : std::ofstream::out);
    if (!output)
        die("Cannot open file " + filename + " to write");
    results_to_wolfram(context, result, output);
    output.close();
    std::cout << "Stored to " + filename << std::endl << std::endl;

    return EXIT_SUCCESS;
}
