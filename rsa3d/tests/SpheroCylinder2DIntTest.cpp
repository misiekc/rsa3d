//
// parallelogramPKua on 08.02.18.
//

#include <fstream>
#include "utility/ShapePairFactory.h"
#include "SpheroCylinder2DIntTest.h"
#include "../Config.h"
#include "../shapes/SpheroCylinder2D.h"
#include "utility/MockBC.h"
#include "../ShapeFactory.h"


namespace
{
    void die(const std::string & _reason)
    {
        std::cerr << _reason << std::endl;
        exit(EXIT_FAILURE);
    }
}

void spheroc_inttest::Context::load(const std::string &_filename) {
    std::ifstream input(_filename);
    if (!input)
        die("Error opening " + std::string(_filename) + " file to read");
    auto config = Config::parse(input);
    this->firstPos = Vector<2>{{config->getDouble("first_x"), config->getDouble("first_y")}};
    this->firstAngle = config->getDouble("first_angle") * M_PI / 180;
    this->secondAngle = config->getDouble("second_angle") * M_PI / 180;
    this->ratio = config->getDouble("ratio");
    this->box_width = config->getDouble("box_width");
    this->box_height = config->getDouble("box_height");
    this->points = config->getDouble("points");
    this->append = config->getString("append") == "true";
    delete config;
}

void spheroc_inttest::perform(const Context &_context, std::vector<spheroc_inttest::Point> &_result) {
    RND rnd;
    MockBC mockBC;

    SpheroCylinder2D first(_context.firstAngle);
    first.vectorTranslate(_context.firstPos);
    SpheroCylinder2D second(_context.secondAngle);

    std::cout << "Generating " << _context.points << " points from " << _context.box_width << " x ";
    std::cout << _context.box_height << " box for:" << std::endl;
    std::cout << first.toString() << std::endl << second.toString() << std::endl;
    for (std::size_t i = 0; i < _context.points; i++) {
        Point point;
        Vector<2> secondTranslation{{
            (rnd.nextValue() - 0.5) * _context.box_width,
            (rnd.nextValue() - 0.5) * _context.box_height}};
        point.position = _context.firstPos + secondTranslation;

        second = SpheroCylinder2D(_context.secondAngle);
        second.vectorTranslate(point.position);
        if (first.overlap(&mockBC, &second))
            point.color = "Red";
        else
            point.color = "Green";
        _result.push_back(point);
    }
}

void spheroc_inttest::results_to_wolfram(const Context &_context, const std::vector<Point> &_results,
                                         std::ostream &_out) {
    // Find Red point with max X coordintate to place the second SC2D
    auto comparator = [](const Point & p1, const Point & p2) {
        if (p1.color == "Green" && p2.color == "Red")
            return true;
        else if (p1.color == "Red" && p2.color == "Green")
            return false;
        else
            return p1.position[0] < p2.position[0];
    };
    Point maxX = *std::max_element(_results.begin(), _results.end(), comparator);

    SpheroCylinder2D first(_context.firstAngle);
    first.vectorTranslate(_context.firstPos);
    SpheroCylinder2D second(_context.secondAngle);
    second.vectorTranslate(maxX.position);

    _out << std::fixed;
    _out << "points = {" << std::endl;
    for (Point point : std::vector<Point>(_results.begin(), _results.end() - 1))
        _out << "    {" << point.color << ", Point[" << point.position << "]}," << std::endl;
    _out << "    {" << _results.back().color << ", Point[" << _results.back().position << "]}};" << std::endl;
    _out << "sc1 = " << first.toWolfram() << ";" << std::endl;
    _out << "sc2 = " << second.toWolfram() << ";" << std::endl;
    _out << std::endl;
    _out << "Graphics[Join[{{Blue, sc2}, {Black, sc1}}, points]]" << std::endl << std::endl;
}

int spheroc_inttest::main(int argc, char **argv) {
    if (argc < 3)
        die("Usage: ./rsa spheroc_inttest [input] [output = spheroc_inttest.nb]");

    Context context;
    context.load(argv[2]);
    ShapeFactory::initShapeClass("SpheroCylinder2D", std::to_string(context.ratio));
    std::vector<Point> results;
    perform(context, results);

    std::string filename = (argc == 4 ? argv[3] : "spheroc_inttest.nb");
    std::ofstream output(filename, context.append ? std::ofstream::app : std::ofstream::out);
    if (!output)
        die("Cannot open file " + filename + " to write");
    results_to_wolfram(context, results, output);
    output.close();
    std::cout << "Stored to " + filename << std::endl << std::endl;

    return EXIT_SUCCESS;
}
