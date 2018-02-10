//
// parallelogramPKua on 08.02.18.
//

#include <fstream>
#include "utility/ShapePairFactory.h"
#include "SpheroCylinder2DExclusionTest.h"
#include "../Config.h"
#include "../shapes/SpheroCylinder2D.h"
#include "utility/MockBC.h"
#include "../ShapeFactory.h"

static const char *const FIRST_SHAPE_COLOR = "Black";
static const char *const SECOND_SHAPE_COLOR = "Blue";
static const char *const INSIDE_COLOR = "Red";
static const char *const OUTSIDE_COLOR = "Green";
static const char *const CIRCLE_FIX_COLOR = "Purple";
static const char *const FIRST_ZONE_COLOR = "Orange";
static const char *const SECOND_ZONE_COLOR = "Cyan";

namespace
{
    void die(const std::string & _reason)
    {
        std::cerr << _reason << std::endl;
        exit(EXIT_FAILURE);
    }
}

void spheroc_extest::Context::load(const std::string &_filename) {
    std::ifstream input(_filename);
    if (!input)
        die("Error opening " + std::string(_filename) + " file to read");
    auto config = Config::parse(input);
    this->firstPos = Vector<2>{{config->getDouble("first_x"), config->getDouble("first_y")}};
    this->firstAngle = config->getDouble("first_angle") * M_PI / 180;
    this->fromAngle = config->getDouble("from_angle") * M_PI / 180;
    this->toAngle = config->getDouble("to_angle") * M_PI / 180;
    this->ratio = config->getDouble("ratio");
    this->box_width = config->getDouble("box_width");
    this->box_height = config->getDouble("box_height");
    this->points = config->getDouble("points");
    this->append = config->getString("append") == "true";
    this->drawSecondShaped = config->getString("draw_second") == "true";
    delete config;
}

void spheroc_extest::perform_inttest(const Context &_context, std::vector<spheroc_extest::Point> &_result) {
    RND rnd;
    MockBC mockBC;

    SpheroCylinder2D first(_context.firstAngle);
    first.vectorTranslate(_context.firstPos);
    SpheroCylinder2D second(_context.fromAngle);

    std::cout << "[OVERLAP TEST] Generating " << _context.points << " points from " << _context.box_width << " x ";
    std::cout << _context.box_height << " box for:" << std::endl;
    std::cout << first.toString() << std::endl << second.toString() << std::endl;

    for (std::size_t i = 0; i < _context.points; i++) {
        Point point = _context.randomPoint(&rnd);
        second = SpheroCylinder2D(_context.fromAngle);
        second.vectorTranslate(point.position);

        if (first.overlap(&mockBC, &second))
            point.color = INSIDE_COLOR;
        else
            point.color = OUTSIDE_COLOR;
        _result.push_back(point);
    }
}

void spheroc_extest::perform_pitest(const Context &_context, std::vector<spheroc_extest::Point> &_result) {
    MockBC mockBC;
    RND rnd;

    SpheroCylinder2D first(_context.firstAngle);
    first.vectorTranslate(_context.firstPos);
    SpheroCylinder2D second(_context.fromAngle);

    std::cout << "[POINTINSIDE TEST] Generating " << _context.points << " points from " << _context.box_width << " x ";
    std::cout << _context.box_height << " box for:" << std::endl;
    std::cout << first.toString() << std::endl << second.toString() << std::endl;

    for (std::size_t i = 0; i < _context.points; i++) {
        Point point = _context.randomPoint(&rnd);
        double pointPosArray[2];
        point.position.copyToArray(pointPosArray);

        if (first.pointInside(&mockBC, pointPosArray, _context.fromAngle, _context.toAngle)) {
            point.color = INSIDE_COLOR;
        } else {
            SpheroCylinder2D angleFromSC(_context.fromAngle);
            angleFromSC.vectorTranslate(point.position);
            SpheroCylinder2D angleToSC(_context.toAngle);
            angleToSC.vectorTranslate(point.position);
            auto withinAngleFromZone = (bool)first.overlap(&mockBC, &angleFromSC);
            auto withinAngleToZone = (bool)first.overlap(&mockBC, &angleToSC);

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

spheroc_extest::Point spheroc_extest::Context::randomPoint(RND * rnd) const {
    Point point;
    Vector<2> secondTranslation{{
            (rnd->nextValue() - 0.5) * box_width,
            (rnd->nextValue() - 0.5) * box_height}};
    point.position = firstPos + secondTranslation;
    return point;
}

void spheroc_extest::results_to_wolfram(const Context &_context, const std::vector<Point> &_results,
                                         std::ostream &_out) {
    SpheroCylinder2D first(_context.firstAngle);
    first.vectorTranslate(_context.firstPos);
    SpheroCylinder2D second(_context.fromAngle);
    second.vectorTranslate(findMaxXRedPoint(_results).position);

    _out << std::fixed;
    _out << "points = {" << std::endl;
    for (Point point : std::vector<Point>(_results.begin(), _results.end() - 1))
        _out << "    {" << point.color << ", Point[" << point.position << "]}," << std::endl;
    _out << "    {" << _results.back().color << ", Point[" << _results.back().position << "]}};" << std::endl;
    _out << "sc1 = " << first.toWolfram() << ";" << std::endl;
    if (_context.drawSecondShaped) {
        _out << "sc2 = " << second.toWolfram() << ";" << std::endl;
        _out << std::endl;
        _out << "Graphics[Join[{{" << SECOND_SHAPE_COLOR << ", sc2}, {" << FIRST_SHAPE_COLOR << ", sc1}}, points]]";
    } else {
        _out << std::endl;
        _out << "Graphics[Join[{{" << FIRST_SHAPE_COLOR << ", sc1}}, points]]";
    }
    _out << std::endl << std::endl;
}

spheroc_extest::Point spheroc_extest::findMaxXRedPoint(const std::vector<spheroc_extest::Point> &_results) {
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

int spheroc_extest::main(int argc, char **argv) {
    if (argc < 3)
        die("Usage: ./rsa spheroc_extest [input] [output = spheroc_extest.nb]");

    Context context;
    context.load(argv[2]);
    std::vector<Point> results;

    std::string mode(argv[1]);
    if (mode == "spheroc_inttest") {
        ShapeFactory::initShapeClass("SpheroCylinder2D", std::to_string(context.ratio));
        perform_inttest(context, results);
    } else if (mode == "spheroc_pitest") {
        ShapeFactory::initShapeClass("SpheroCylinder2D", std::to_string(context.ratio));
        perform_pitest(context, results);
    }

    std::string filename = (argc == 4 ? argv[3] : "spheroc_inttest.nb");
    std::ofstream output(filename, context.append ? std::ofstream::app : std::ofstream::out);
    if (!output)
        die("Cannot open file " + filename + " to write");
    results_to_wolfram(context, results, output);
    output.close();
    std::cout << "Stored to " + filename << std::endl << std::endl;

    return EXIT_SUCCESS;
}
