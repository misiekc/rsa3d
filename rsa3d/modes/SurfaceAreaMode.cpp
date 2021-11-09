//
// Created by pkua on 06.11.2021.
//

#include <sstream>

#include "SurfaceAreaMode.h"
#include "../PackingGenerator.h"

void SurfaceAreaMode::run() {
    PackingGenerator pg(0, 0, &this->params);
    double surfaceArea = pg.getSurface().getRealArea();
    std::stringstream out;
    out.precision(std::numeric_limits<double>::max_digits10);
    out << surfaceArea;
    std::cout << out.rdbuf() << std::endl;
}

void SurfaceAreaMode::printHelp(std::ostream &out, const ProgramArguments &arguments) {
    out << arguments.formatUsage("") << std::endl;
    out << std::endl;
    out << "For given input file, it calculates the collector surface area (or its d-dimensional" << std::endl;
    out << "counterpart). Note, that this value is non-trivial for example for curved surfaces." << std::endl;
}
