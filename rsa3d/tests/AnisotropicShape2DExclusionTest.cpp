//
// PKua on 08.02.18.
//

#include <fstream>
#include <bits/unique_ptr.h>
#include <memory>
#include "utility/ShapePairFactory.h"
#include "AnisotropicShape2DExclusionTest.h"
#include "../Config.h"
#include "../shapes/SpheroCylinder2D.h"
#include "utility/MockBC.h"
#include "../ShapeFactory.h"
#include "../Utils.h"

static const char *const FIRST_SHAPE_COLOR = "Black";
static const char *const SECOND_SHAPE_COLOR = "Blue";
static const char *const INSIDE_COLOR = "Red";
static const char *const OUTSIDE_COLOR = "Green";
static const char *const CIRCLE_FIX_COLOR = "Purple";
static const char *const FIRST_ZONE_COLOR = "Orange";
static const char *const SECOND_ZONE_COLOR = "Cyan";


void as2d_extest::Context::load(const std::string &_filename) {
    std::ifstream input(_filename);
    if (!input)
        die("Error opening " + std::string(_filename) + " file to read");
    auto config = Config::parse(input);
    this->particle = config->getString("particle");
    this->attributes = config->getString("attr");
    this->firstPos = Vector<2>{{config->getDouble("first_x"), config->getDouble("first_y")}};
    this->firstAngle = config->getDouble("first_angle") * M_PI / 180;

    std::string modeStr = config->getString("mode");
    if (modeStr == "overlap")
        this->mode = Mode::OVERLAP;
    else if (modeStr == "point_inside")
        this->mode = Mode::POINT_INSIDE;
    else
        throw std::runtime_error("Unknown test mode: " + modeStr);

    if (this->mode == Mode::OVERLAP) {
        this->secondAngle = config->getDouble("second_angle") * M_PI / 180;
    } else if (this->mode == Mode::POINT_INSIDE) { // Mode::POINT_INSIDE
        this->fromAngle = config->getDouble("from_angle") * M_PI / 180;
        this->toAngle = config->getDouble("to_angle") * M_PI / 180;
    }

    this->box_width = config->getDouble("box_width");
    this->box_height = config->getDouble("box_height");
    this->points = config->getDouble("points");
    this->append = config->getString("append") == "true";
    delete config;

    std::cout << "Context loaded from " << _filename << std::endl;
}

void as2d_extest::perform_inttest(const Context &_context, std::vector<as2d_extest::Point> &_result) {
    RND rnd;
    MockBC mockBC;

    _context.printInfo();
    auto first = generate_shape(_context.firstAngle, _context.firstPos);
    for (std::size_t i = 0; i < _context.points; i++) {
        Point point = _context.randomPoint(&rnd);
        auto second = generate_shape(_context.secondAngle, point.position);

        if (first->overlap(&mockBC, second.get()))
            point.color = INSIDE_COLOR;
        else
            point.color = OUTSIDE_COLOR;
        _result.push_back(point);
    }
}

std::unique_ptr<AnisotropicShape2D> as2d_extest::generate_shape(double angle, const Vector<2> &position) {
    auto shape = (AnisotropicShape2D*)ShapeFactory::createShape(nullptr);
    double posAngle[3] = {position[0], position[1], angle};
    shape->translate(posAngle);
    shape->rotate(posAngle + 2);
    return std::unique_ptr<AnisotropicShape2D>(shape);
}

void as2d_extest::perform_pitest(const Context &_context, std::vector<as2d_extest::Point> &_result) {
    MockBC mockBC;
    RND rnd;

    _context.printInfo();
    auto first = generate_shape(_context.firstAngle, _context.firstPos);
    for (std::size_t i = 0; i < _context.points; i++) {
        Point point = _context.randomPoint(&rnd);
        double pointPosArray[2];
        point.position.copyToArray(pointPosArray);

        if (first->pointInside(&mockBC, pointPosArray, _context.fromAngle, _context.toAngle)) {
            point.color = INSIDE_COLOR;
        } else {
            auto angleFromShape = generate_shape(_context.fromAngle, point.position);
            auto angleToShape = generate_shape(_context.toAngle, point.position);
            auto withinAngleFromZone = (bool)first->overlap(&mockBC, angleFromShape.get());
            auto withinAngleToZone = (bool)first->overlap(&mockBC, angleToShape.get());

            if (withinAngleFromZone && withinAngleToZone)
                point.color = CIRCLE_FIX_COLOR;
            else if (withinAngleFromZone)
                point.color = FIRST_ZONE_COLOR;
            else if (withinAngleToZone)
                point.color = SECOND_ZONE_COLOR;
            else
                point.color = OUTSIDE_COLOR;
        }
        _result.push_back(point);
    }
}

as2d_extest::Point as2d_extest::Context::randomPoint(RND * rnd) const {
    Point point;
    Vector<2> secondTranslation{{
            (rnd->nextValue() - 0.5) * box_width,
            (rnd->nextValue() - 0.5) * box_height}};
    point.position = firstPos + secondTranslation;
    return point;
}

void as2d_extest::Context::printInfo() const {
    if (this->mode == Mode::OVERLAP)
        std::cout << "[OVERLAP TEST] ";
    else    // Mode::POINT_INSIDE
        std::cout << "[POINT INSIDE TEST] ";

    auto first = generate_shape(this->firstAngle, this->firstPos);
    std::cout << "Generating " << this->points << " points from " << this->box_width << " x ";
    std::cout << this->box_height << " box for:" << std::endl;
    std::cout << first->toString() << std::endl;

    if (this->mode == Mode::OVERLAP)
        std::cout << "Second angle: " << this->secondAngle * 180 / M_PI << std::endl;
    else    // Mode::POINT_INSIDE
        std::cout << "Angle range: [" << this->fromAngle * 180 / M_PI << ", " << this->toAngle * 180 / M_PI << "]" << std::endl;
}

void as2d_extest::results_to_wolfram(const Context &_context, const std::vector<Point> &_results,
                                         std::ostream &_out) {
    if (_context.append)
        _out << std::endl << std::endl;

    _out << std::fixed;
    _out << "points = {" << std::endl;
    for (Point point : std::vector<Point>(_results.begin(), _results.end() - 1))
        _out << "    {" << point.color << ", Point[" << point.position << "]}," << std::endl;
    _out << "    {" << _results.back().color << ", Point[" << _results.back().position << "]}};" << std::endl;

    auto first = generate_shape(_context.firstAngle, _context.firstPos);
    _out << "shape1 = " << first->toWolfram() << ";" << std::endl;
    if (_context.mode == Context::Mode::OVERLAP) {
       auto second = generate_shape(_context.secondAngle, find_max_x_red_point(_results).position);
        _out << "shape2 = " << second->toWolfram() << ";" << std::endl << std::endl;
        _out << "Graphics[Join[{{" << SECOND_SHAPE_COLOR << ", shape2}, {" << FIRST_SHAPE_COLOR << ", shape1}}, points]]";
    } else {  // Context::Mode::POINT_INSIDE
        _out << std::endl << "Graphics[Join[{{" << FIRST_SHAPE_COLOR << ", shape1}}, points]]";
    }
}

as2d_extest::Point as2d_extest::find_max_x_red_point(const std::vector<as2d_extest::Point> &_results) {
    auto comparator = [](const Point & p1, const Point & p2) {
        if (p1.color != INSIDE_COLOR && p2.color == INSIDE_COLOR)
            return true;
        else if (p1.color == INSIDE_COLOR && p2.color != INSIDE_COLOR)
            return false;
        else
            return p1.position[0] < p2.position[0];
    };

    return *max_element(_results.begin(), _results.end(), comparator);
}

int as2d_extest::main(int argc, char **argv) {
    if (argc < 3)
        die("Usage: ./rsa as2d_extest [input] [output = out.nb]");

    Context context;
    context.load(argv[2]);
    std::vector<Point> results;

    ShapeFactory::initShapeClass(context.particle, context.attributes);
    if (context.mode == Context::Mode::OVERLAP)
        perform_inttest(context, results);
    else  // Context::Mode::POINT_INSIDE
        perform_pitest(context, results);

    std::string filename = (argc == 4 ? argv[3] : "out.nb");
    std::ofstream output(filename, context.append ? std::ofstream::app : std::ofstream::out);
    if (!output)
        die("Cannot open file " + filename + " to write");
    results_to_wolfram(context, results, output);
    output.close();
    std::cout << "Stored to " + filename << std::endl << std::endl;

    return EXIT_SUCCESS;
}
