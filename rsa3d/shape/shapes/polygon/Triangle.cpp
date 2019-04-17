//
// Created by pkua on 10.04.19.
//

#include "Triangle.h"

void Triangle::initClass(const std::string &args) {
    std::string attr = Triangle::preparePolygonAttributes(args);
    Polygon::initClass(attr);
}

/* Separate function so it is testable */
std::string Triangle::preparePolygonAttributes(const std::string &triangleAttr) {
    double a, b, c;
    Triangle::parseAttributes(triangleAttr, a, b, c);
    Vector<2> thirdVertex = Triangle::calculateThirdVertex(a, b, c);

    std::ostringstream out;
    out << "3 xy " << -c/2 << " 0 " << c/2 << " 0 " << thirdVertex[0] << " " << thirdVertex[1];
    out << " 3 0 1 1 2 2 0 starHelperSegments";
    return out.str();
}

void Triangle::parseAttributes(const std::string &args, double &a, double &b, double &c) {
    std::istringstream in(args);
    in >> a >> b >> c;
    if (!in) throw std::domain_error("Expected 3 doubles as side lengths");

    // Ensure order a <= b <= c
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > c) std::swap(a, c);

    Validate(a > 0);
    Validate(b > 0);
    Validate(c > 0);
    Validate(a + b > c);        // Triangle inequality
}

/* Calculates the coordinates of the third triangle vertex assuming that its base is on X axis. It assumes that
 * positions of two first vertices are (-c/2, 0), (c/2, 0) and the third one has y > 0. Parameters:
 * a - the length of the right side
 * b - the length of the left side
 * c - the length of rhe horizontal side lying on x axis
 */
Vector<2> Triangle::calculateThirdVertex(double a, double b, double c) {
    double x = (c*c + a*a - b*b)/(2*c);
    return {{c/2 - x, std::sqrt(a*a - x*x)}};
}