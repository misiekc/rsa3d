//
// Created by Kasperek, Wojciech on 10/03/2018.
//

#include "Rectangle.h"


double Rectangle::getXOnCircle(double cx, double cy, double r, double angle) {
    return cx + cos(angle) * r;
}

double Rectangle::getYOnCircle(double cx, double cy, double r, double angle) {
    return cy + sin(angle) * r;
}

void Rectangle::setAngle(double angle) {
    this->orientation[0] = angle;
}

double Rectangle::getVoxelSize() {
    return Rectangle::voxelSize;
}

double Rectangle::getNeighbourListCellSize() {
    return Rectangle::neighbourListCellSize;
}

void Rectangle::initClass(const std::string &args) {
    double ratio = std::stod(args);
    double s = sqrt(1.0/ratio);
    if (ratio >= 0) {
        Rectangle::longer = s * ratio;
        Rectangle::shorter = s;
    } else {
        Rectangle::longer = s;
        Rectangle::shorter = s * ratio;
    }
    Rectangle::voxelSize = Rectangle::shorter / 2;
    Rectangle::neighbourListCellSize = Rectangle::longer * 2;
}

double Rectangle::getVolume() {
    return 1;
}

Shape<2, 1> *Rectangle::create(RND *rnd) {
    return new Rectangle();
}

Rectangle::Rectangle() {
    this->a = Rectangle::longer;
    this->b = Rectangle::shorter;
    this->halfA = a/2;
    this->halfB = b/2;
    this->setAngle(0);
}

Rectangle::Rectangle(const Rectangle &other) {
    this->a = other.a;
    this->b = other.b;
    this->halfA = other.halfA;
    this->halfB = other.halfB;
    this->setAngle(other.getAngle());
}

Rectangle &Rectangle::operator=(const Rectangle &el) {
    this->a = el.a;
    this->b = el.b;
    this->halfA = el.halfA;
    this->halfB = el.halfB;
    this->setAngle(el.getAngle());
    return *this;
}

void Rectangle::rotate(double *v) {
    Shape::rotate(v);
}

int Rectangle::pointInside(BoundaryConditions *bc, double *da) {
    // apply bc
    double ta[2];
    double tmp[2];

    tmp[0] = da[0]; tmp[1] = da[1];
    bc->getTranslation(ta, this->position, tmp);
    da[0] += ta[0];
    da[1] += ta[1];

    // rotate point relatively to rectangle center
    Vector<2> pointRotated = rotatePoint(Vector<2>(da), getAntiRotationMatrix(), Vector<2>(position));
    if (!isInsideRect(pointRotated, position[0], position[1], a + b, b + b)) {
        return 0;
    }

    // point is no further than halfB from rect's side
    // but can be further than halfB from rect's corner

    double cx = position[0] + halfA;
    double cy = position[1] + halfB;
    double r = halfB;
    if (pointRotated[0] > cx && pointRotated[1] > cy) { // corner 1,1
        return isInsideCircle(pointRotated, cx, cy, r);
    }

    cx = position[0] + halfA;
    cy = position[1] - halfB;
    if (pointRotated[0] > cx && pointRotated[1] < cy) { // corner 1,-1
        return isInsideCircle(pointRotated, cx, cy, r);
    }

    cx = position[0] - halfA;
    cy = position[1] - halfB;
    if (pointRotated[0] < cx && pointRotated[1] < cy) { // corner -1,-1
        return isInsideCircle(pointRotated, cx, cy, r);
    }

    cx = position[0] - halfA;
    cy = position[1] + halfB;
    if (pointRotated[0] < cx && pointRotated[1] > cy) { // corner -1,-1
        return isInsideCircle(pointRotated, cx, cy, r);
    }

    // not in a corner, but in a excluding rect
    return 1;
}


int Rectangle::pointInside(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) {
    // let's have it ordered
    if (angleFrom > angleTo) {
        double tmp = angleFrom;
        angleFrom = angleTo;
        angleTo = angleFrom;
    }

    // and normalized
    normalizeAngleRange(&angleFrom, &angleTo, M_PI);

    // 90 degrees means minimal exclusion zone
    if (angleTo - angleFrom > M_PI) {
        return pointInside(bc, da);
    }

    double pi2 = M_PI / 2;
    // angle range is M_PI max
    int countFrom = (int) floor(angleFrom / pi2);
    int countTo = (int) floor(angleTo / pi2);
    if (countTo - countFrom == 0) { // cool, easy
        // let's move on
    } else if (countTo - countFrom == 1) { // double check
        return pointInside(bc, da, angleFrom, pi2 * countTo) && pointInside(bc, da, pi2 * countTo, angleTo);
    } else { // triple check
        return pointInside(bc, da, angleFrom, pi2 * (countFrom + 1)) && pointInside(bc, da, pi2 * (countFrom + 1), pi2 * countTo)
               && pointInside(bc, da, pi2 * countTo, angleTo);
    }

    // apply bc
    double ta[2];
    double tmp[2];

    tmp[0] = da[0]; tmp[1] = da[1];
    bc->getTranslation(ta, this->position, tmp);
    da[0] += ta[0];
    da[1] += ta[1];

    // rotate point relatively to rectangle center
    Vector<2> pointRotated = rotatePoint(Vector<2>(da), getAntiRotationMatrix(), Vector<2>(position));
    double x = pointRotated[0];
    double y = pointRotated[1];

    // if it is not in excluding rectangle, we can quit
    if (!checkExcludingRectangle(angleFrom, angleTo, x, y)) {
        return 0;
    }

    // there are two options:
    bool option = countFrom % 2 == 1;

    // 1, 1
    // setup
    double cx = halfA;
    double cy = halfB;
    double r = option ? b : a;
    double angle1 = angleFrom - option * 90;
    double angle2 = angleTo - option * 90;
    if (!checkTangents(x, y, cx, cy, r, angle1, angle2)) {
        return false;
    }

    // -1, 1
    cx = -halfA;
    cy = halfB;
    r = option ? a : b;
    angle1 = angleFrom + 90 - option * 90;
    angle2 = angleTo + 90 - option * 90;
    if (!checkTangents(x, y, cx, cy, r, angle1, angle2)) {
        return false;
    }

    // -1, -1
    cx = -halfA;
    cy = -halfB;
    r = option ? b : a;
    angle1 = angleFrom + 180 - option * 90;
    angle2 = angleTo + 180 - option * 90;
    if (!checkTangents(x, y, cx, cy, r, angle1, angle2)) {
        return false;
    }

    // 1, -1
    cx = halfA;
    cy = -halfB;
    r = option ? a : b;
    angle1 = angleFrom + 270 - option * 90;
    angle2 = angleTo + 270 - option * 90;
    if (!checkTangents(x, y, cx, cy, r, angle1, angle2)) {
        return false;
    }

    return true;
}

bool Rectangle::checkTangents(double x, double y, double cx, double cy, double r, double angle1, double angle2) {// calculate tangents
    double p1x = getXOnCircle(cx, cy, r, angle1);
    double p1y = getYOnCircle(cx, cy, r, angle2);

    double a1 = p1x - cx;
    double b1 = p1y - cy;
    double c1 = a1 * (-cx) + b1 * (-cy) - r * r;

    double p2x = getXOnCircle(cx, cy, r, angle1);
    double p2y = getYOnCircle(cx, cy, r, angle2);

    double a2 = p2x - cx;
    double b2 = p2y - cy;
    double c2 = a2 * (-cx) + b2 * (-cy) - r * r;

    double f1 = a1 * x + b1 * y + c1;
    double f2 = a2 * x + b2 * y + c2;

    // check if "inside" tangents
    bool isIn = f1 < 0 && f2 < 0;

    // check circle
    if (isIn && x >= std::__1::min(p2x, p1x) && x <= std::__1::max(p2x, p1x) && y >= std::__1::min(p2y, p1y) && y <= std::__1::max(p2y, p1y)) {
        isIn = (cx - x) * (cx - x) + (cy - y) * (cy - y) < r * r;
    }
    return isIn;
}

bool Rectangle::checkExcludingRectangle(double angleFrom, double angleTo, double x, double y) const {// so we have angle range in 0-90, 90-180 or whatever else
    // let's first calculate "excluding rectangle"
    Rectangle rec = Rectangle(*this);
    rec.setAngle(angleFrom);
    Matrix<2, 2> rotationFrom = rec.getRotationMatrix();
    rec.setAngle(angleTo);
    Matrix<2, 2> rotationTo = rec.getRotationMatrix();
    double p[2] = {halfA, halfB};
    Vector<2> newSizeFrom = Vector<2>();
    Vector<2> newSizeTo = Vector<2>();
    getNewSize(rotationFrom, p, newSizeFrom);
    p[0] = halfA;
    p[1] = halfB;
    getNewSize(rotationTo, p, newSizeTo);

    double newA = (newSizeFrom[0] > newSizeTo[0] ? newSizeFrom[0] : newSizeTo[0]) + a;
    double newB = (newSizeFrom[1] > newSizeTo[1] ? newSizeFrom[1] : newSizeTo[1]) + b;


    // if the point is outside rectangle, all work is done
    if (x >= newA || x <= -newA || y >= newB || y <= -newB) {
        return 0;
    }
    return 1;
}

void Rectangle::getNewSize(const Matrix<2, 2, double> &rotationFrom, const double *p, const Vector<2, double> &newSizeFrom) const {
    Vector<2> helpMe = Vector<2>(p);
    Vector<2> result = rotationFrom * helpMe;
    newSizeFrom[0] = result[0];
    newSizeFrom[1] = result[1];
    helpMe[0] -= a;
    result = rotationFrom * helpMe;
    if (result[0] > newSizeFrom[0]) {
        newSizeFrom[0] = result[0];
    }
    if (result[1] > newSizeFrom[1]) {
        newSizeFrom[1] = result[1];
    }
    helpMe[1] -= b;
    result = rotationFrom * helpMe;
    if (result[0] > newSizeFrom[0]) {
        newSizeFrom[0] = result[0];
    }
    if (result[1] > newSizeFrom[1]) {
        newSizeFrom[1] = result[1];
    }
    helpMe[0] += a;
    result = rotationFrom * helpMe;
    if (result[0] > newSizeFrom[0]) {
        newSizeFrom[0] = result[0];
    }
    if (result[1] > newSizeFrom[1]) {
        newSizeFrom[1] = result[1];
    }
}


int Rectangle::isInsideCircle(const Vector<2, double> &point, double cx, double cy, double r) const {
    if ((cx - point[0]) * (cx - point[0]) + (cy - point[1]) * (cy - point[1]) <= r * r) {
        return 1;
    }
    return 0;
}

int Rectangle::isInsideRect(const Vector<2, double> &point, double cx, double cy, double a, double b) {
    if (abs(cx - point[0]) <= a/2 && abs(cy - point[1]) <= b/2) {
        return 1;
    }
    return 0;
}

Vector<2> Rectangle::rotatePoint(const Vector<2>& point, const Matrix<2,2>& rotation) {
    return rotation * point;
}

Vector<2> Rectangle::rotatePoint(const Vector<2>& point, const Matrix<2,2>& rotation, const Vector<2>& center) {
    Vector<2> translatedPoint = Vector();
    translatedPoint[0] = point[0] - center[0];
    translatedPoint[1] = point[1] - center[1];
    Vector<2> rotated =  rotation * translatedPoint;
    rotated[0] += center[0];
    rotated[1] += center[1];
    return rotated;
}

int Rectangle::overlap(BoundaryConditions *bc, Shape<2, 1> *s) {
    Rectangle rectangle = *((Rectangle *)s);
    this->applyBC(bc, &rectangle);

    // let's say our rectangle stays aligned to axis
    // rotate the other one
    // and check if corners are inside
    Matrix<2,2> rotation = getAntiRotationMatrix() * rectangle.getRotationMatrix();
    Vector<2> rectPosition = Vector<2>(rectangle.position);
    double x = position[0];
    double y = position[1];
    double p[2] = {rectangle.position[0] + halfA, rectangle.position[1] + halfB}; // 1,1
    Vector<2> point = rotatePoint(Vector<2>(p), rotation, rectPosition);
    if (isInsideRect(point, x, y, a, b)) {
        return 1;
    }

    p[0] -= a; // -1, 1
    point = rotatePoint(Vector<2>(p), rotation, rectPosition);
    if (isInsideRect(point, x, y, a, b)) {
        return 1;
    }

    p[1] -= b; // -1, -1
    point = rotatePoint(Vector<2>(p), rotation, rectPosition);
    if (isInsideRect(point, x, y, a, b)) {
        return 1;
    }

    p[0] += a; // 1, -1
    point = rotatePoint(Vector<2>(p), rotation, rectPosition);
    if (isInsideRect(point, x, y, a, b)) {
        return 1;
    }

    // okay, now other way round
    // second rectangle stays aligned, our is rotated
    rotation = rectangle.getAntiRotationMatrix() * getRotationMatrix();

    rectPosition = Vector<2>(position);
    x = rectangle.position[0];
    y = rectangle.position[1];
    p[0] = halfA + position[0];
    p[1] = halfB + position[1]; // 1,1
    point = rotatePoint(Vector<2>(p), rotation, rectPosition);
    if (isInsideRect(point, x, y, a, b)) {
        return 1;
    }

    p[0] -= a; // -1, 1
    point = rotatePoint(Vector<2>(p), rotation, rectPosition);
    if (isInsideRect(point, x, y, a, b)) {
        return 1;
    }

    p[1] -= b; // -1, -1
    point = rotatePoint(Vector<2>(p), rotation, rectPosition);
    if (isInsideRect(point, x, y, a, b)) {
        return 1;
    }

    p[0] += a; // 1, -1
    point = rotatePoint(Vector<2>(p), rotation, rectPosition);
    if (isInsideRect(point, x, y, a, b)) {
        return 1;
    }

    // no overlap
    return 0;
}

std::string Rectangle::toWolfram() const {
    std::stringstream out;

    Vector<2> thisPos(this->position);
    out << std::fixed;
    out << "GeometricTransformation[Rectangle[{" << -halfA <<"," << -halfB << "}, {" << halfA << "," << halfB << "}]," << std::endl;
    out << "    {RotationMatrix[" << this->getAngle() << "], " << thisPos << "}]";

    return out.str();
}


std::string Rectangle::toString() {
    std::stringstream out;

    Vector<2> thisPos(this->position);
    out << "Rectangle{position: " << thisPos << "; a: " << this->a;
    out << "; b: " << this->b << "; angle: " << this->getAngle() << "}";
    return out.str();
}

std::string Rectangle::toPovray() const{
    // ignore
    return "";
}

void Rectangle::store(std::ostream &f) const{
    Shape::store(f);
    f.write((char *)(&this->a), sizeof(double));
    f.write((char *)(&this->b), sizeof(double));
}

void Rectangle::restore(std::istream &f){
    Shape::restore(f);
    f.read((char *)(&this->a), sizeof(double));
    f.read((char *)(&this->b), sizeof(double));
    this->halfA = a/2;
    this->halfB = b/2;
}