//
// Created by pkua on 10.04.19.
//

#include "Triangle.h"

void Triangle::initClass(const std::string &args) {
    double a, b, c;
    std::istringstream in(args);
    in >> a >> b >> c;
    Validate(a > 0);
    Validate(b > 0);
    Validate(c > 0);

    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    Validate(a + b > c);        // Triangle inequality

    Vector<2> thirdVertex = Triangle::calculateThirdVertex(a, b, c);

    std::ostringstream out;
    out << "3 xy -0.5 0 0.5 0 " << thirdVertex[0] << " " << thirdVertex[1] << " 3 0 1 1 2 2 0 starHelperSegments";
    std::cout << out.str() << std::endl;
    Polygon::initClass(out.str());
}

Vector<2> Triangle::calculateThirdVertex(double a, double b, double c) {
    double x = (c*c + a*a - b*b)/(2*c);
    return {{c/2 - x, std::sqrt(a*a - x*x)}};
}
