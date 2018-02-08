//
// Created by PKua on 08.02.18.
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
    delete config;
}

void spheroc_inttest::perform(const Context &_context, std::vector<Point> &_result) {
    RND rnd;
    MockBC mockBC;

    SpheroCylinder2D first(_context.firstAngle);
    first.vectorTranslate(_context.firstPos);
    for (std::size_t i = 0; i < _context.points; i++) {
        Point point;
        Vector<2> secondTranslation{{
            (2 * rnd.nextValue() - 1) * _context.box_width,
            (2 * rnd.nextValue() - 1) * _context.box_height}};
        point.position = _context.firstPos + secondTranslation;

        SpheroCylinder2D second(_context.secondAngle);
        second.vectorTranslate(point.position);
        if (first.overlap(&mockBC, &second))
            point.color = "Red";
        else
            point.color = "Green";
        _result.push_back(point);
    }
}

void spheroc_inttest::print_results(const Context &_context, const std::vector<spheroc_inttest::Point> &_results,
                                    std::ostream &_out) {
    SpheroCylinder2D first(_context.firstAngle);
    first.vectorTranslate(_context.firstPos);
    SpheroCylinder2D second(_context.secondAngle);

    _out << "points = {" << std::endl;
    for (Point point : std::vector<Point>(_results.begin(), _results.end() - 1))
        _out << "    {" << point.color << ", Point[" << point.position << "]}," << std::endl;
    _out << "    {" << _results.back().color << ", Point[" << _results.back().position << "]}};" << std::endl;
    _out << "sc1 = " << first.toWolfram() << ";" << std::endl;
    _out << "Graphics[Join[{{Black, sc}}, points]]" << std::endl;
    _out << std::endl;
    _out << "sc2 = " << second.toWolfram() << ";" << std::endl;
    _out << "Graphics[{Black, sc2}]";
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
    std::ofstream output(filename);
    if (!output)
        die("Cannot open file " + filename + " to write");
    print_results(context, results, output);
    output.close();

    return EXIT_SUCCESS;
}
