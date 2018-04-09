//
// Created by Kasperek, Wojciech on 10/03/2018.
//

#include "Rectangle.h"

double Rectangle::neighbourListCellSize;
double Rectangle::voxelSize;
double Rectangle::longer;
double Rectangle::shorter;

Shape<2, 1> *Rectangle::create(RND *rnd) {
    return new Rectangle();
}

Rectangle::Rectangle() : AnisotropicShape2D() {
    this->a = Rectangle::longer;
    this->b = Rectangle::shorter;
    this->halfA = a / 2;
    this->halfB = b / 2;
    this->setAngle(0);
}

double Rectangle::getXOnCircle(double cx, double cy, double r, double angle) const {
    return cx + cos(angle) * r;
}

double Rectangle::getYOnCircle(double cx, double cy, double r, double angle) const {
    return cy + sin(angle) * r;
}

void Rectangle::setAngle(double angle) {
    AnisotropicShape2D::setAngle(angle);

    const double *position = this->getPosition();
    double p[2] = {halfA, halfB}; // 1, 1
    const Matrix<2, 2, double> &rotation = getRotationMatrix();
    Vector<2> point = rotatePoint(Vector<2>(p), rotation);
    xs[0] = point[0] + position[0];
    xs[4] = point[0] + position[0];
    ys[0] = point[1] + position[1];
    ys[4] = point[1] + position[1];
    p[0] -= a; // -1, 1
    point = rotatePoint(Vector<2>(p), rotation);
    xs[1] = point[0] + position[0];
    ys[1] = point[1] + position[1];
    p[1] -= b; // -1, -1
    point = rotatePoint(Vector<2>(p), rotation);
    xs[2] = point[0] + position[0];
    ys[2] = point[1] + position[1];
    p[0] += a; // 1, - 1
    point = rotatePoint(Vector<2>(p), rotation);
    xs[3] = point[0] + position[0];
    ys[3] = point[1] + position[1];
}

double Rectangle::getVoxelSize() const {
    return Rectangle::voxelSize;
}

double Rectangle::getNeighbourListCellSize() const {
    return Rectangle::neighbourListCellSize;
}

void Rectangle::initClass(const std::string &args) {
    double ratio = std::stod(args);
    double s = sqrt(1.0 / ratio);
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

double Rectangle::getVolume() const {
    return 1;
}

double Rectangle::getVoxelAngularSize() const {
    return M_PI;
}

// Shape::translate was made non-virtual and one should override Positioned::setPosition instead (see documentation)
void Rectangle::setPosition(const double *position) {
    Vector<2, double> translation = Vector<2, double>(position) - Vector<2, double>(getPosition());
    for (int i = 0; i < 5; i++) {
        xs[i] += translation[0];
        ys[i] += translation[1];
    }

    Shape::setPosition(position);
}

// never used
int Rectangle::pointInside(BoundaryConditions *bc, double *da) const {
    // apply bc
    double ta[2];
    double tmp[2];
    double *position = this->getPosition();

    tmp[0] = da[0];
    tmp[1] = da[1];
    bc->getTranslation(ta, position, tmp);
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
    if (pointRotated[0] < cx && pointRotated[1] > cy) { // corner -1,1
        return isInsideCircle(pointRotated, cx, cy, r);
    }

    // not in a corner, but in a excluding rect
    return 1;
}


int Rectangle::pointInside(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) const {
    // let's have it ordered
    if (angleFrom > angleTo) {
        double tmp = angleFrom;
        angleFrom = angleTo;
        angleTo = tmp;
    }

    double angle = this->getAngle();
    angleFrom -= angle;
    angleTo -= angle;

    // and normalized
    normalizeAngleRange(&angleFrom, &angleTo, M_PI);

    // 90 degrees means minimal exclusion zone
    if (angleTo - angleFrom >= M_PI) {
        return pointInside(bc, da);
    }

    double pi2 = M_PI / 2;
    // angle range is M_PI max
    int countFrom = (int) floor(angleFrom / pi2);
    int countTo = (int) floor(angleTo / pi2);
    if (countTo - countFrom == 0 || (countTo - countFrom == 1 && angleTo == countTo * pi2)) {
        return pointInsideInternal(bc, da, angleFrom, angleTo);
    } else if (countTo - countFrom == 1) { // double check
        return pointInsideInternal(bc, da, angleFrom, pi2 * countTo)
               && pointInsideInternal(bc, da, pi2 * countTo, angleTo);
    } else { // triple check
        return pointInsideInternal(bc, da, angleFrom, pi2 * (countFrom + 1))
               && pointInsideInternal(bc, da, pi2 * (countFrom + 1), pi2 * countTo)
               && pointInsideInternal(bc, da, pi2 * countTo, angleTo);
    }
}

int Rectangle::pointInsideInternal(BoundaryConditions *bc, double *da, double angleFrom, double angleTo) const {
    double pi2 = M_PI / 2;
    int countFrom = (int) floor(angleFrom / pi2);

    // apply bc
    double ta[2];
    double tmp[2];
    double *position = this->getPosition();

    tmp[0] = da[0];
    tmp[1] = da[1];
    bc->getTranslation(ta, position, tmp);
    Vector<2, double> point = Vector<2>(da);
    point[0] += ta[0];
    point[1] += ta[1];

    // rotate point relatively to rectangle center
    const Matrix<2, 2, double> &rotation = getAntiRotationMatrix();

    Vector<2> pointRotated = rotatePoint(point, rotation, Vector<2>(position));
    double x = pointRotated[0];
    double y = pointRotated[1];

    // if it is not in excluding rectangle, we can quit
    if (!isInsideExcludingRectangle(angleFrom, angleTo, x, y)) {
        return false;
    }

    // inside rect?
    if (isInsideRect(pointRotated, position[0], position[1], a, b)) {
        return true;
    }

    // there are four quarters:
    int quarter = countFrom % 4;

    // 1, 1
    // setup
    double cx = halfA + position[0];
    double cy = halfB + position[1];
    double r = quarter % 2 == 0 ? halfA : halfB;
    double desiredQuarter = 0;
    double angle1 = angleFrom - (quarter - desiredQuarter) * pi2;
    double angle2 = angleTo - (quarter - desiredQuarter) * pi2;
    if (!checkTangents(x, y, cx, cy, r, angle1, angle2)) {
        return false;
    }

    // -1, 1
    cx = -halfA + position[0];
    cy = halfB + position[1];
    r = quarter % 2 == 1 ? halfA : halfB;
    desiredQuarter = 1;
    angle1 = angleFrom - (quarter - desiredQuarter) * pi2;
    angle2 = angleTo - (quarter - desiredQuarter) * pi2;
    if (!checkTangents(x, y, cx, cy, r, angle1, angle2)) {
        return false;
    }

    // -1, -1
    cx = -halfA + position[0];
    cy = -halfB + position[1];
    r = quarter % 2 == 0 ? halfA : halfB;
    desiredQuarter = 2;
    angle1 = angleFrom - (quarter - desiredQuarter) * pi2;
    angle2 = angleTo - (quarter - desiredQuarter) * pi2;
    if (!checkTangents(x, y, cx, cy, r, angle1, angle2)) {
        return false;
    }

    // 1, -1
    cx = halfA + position[0];
    cy = -halfB + position[1];
    r = quarter % 2 == 1 ? halfA : halfB;
    desiredQuarter = 3;
    angle1 = angleFrom - (quarter - desiredQuarter) * pi2;
    angle2 = angleTo - (quarter - desiredQuarter) * pi2;
    if (!checkTangents(x, y, cx, cy, r, angle1, angle2)) {
        return false;
    }

    return true;
}

bool Rectangle::checkTangents(double x, double y, double cx, double cy, double r, double angle1,
                              double angle2) const {
    // calculate tangents
    double p1x = getXOnCircle(cx, cy, r, angle1);
    double p1y = getYOnCircle(cx, cy, r, angle1);

    double a1 = p1x - cx;
    double b1 = p1y - cy;
    double c1 = a1 * (-cx) + b1 * (-cy) - r * r;

    double p2x = getXOnCircle(cx, cy, r, angle2);
    double p2y = getYOnCircle(cx, cy, r, angle2);

    double a2 = p2x - cx;
    double b2 = p2y - cy;
    double c2 = a2 * (-cx) + b2 * (-cy) - r * r;

    double f1 = a1 * x + b1 * y + c1;
    double f2 = a2 * x + b2 * y + c2;

    // check if "inside" tangents
    bool isIn = f1 <= 0 && f2 <= 0;

    // check circle
    if (isIn && x >= std::min(p2x, p1x) && x <= std::max(p2x, p1x) && y >= std::min(p2y, p1y) &&
        y <= std::max(p2y, p1y)) {
        isIn = (cx - x) * (cx - x) + (cy - y) * (cy - y) <= r * r;
    }
    return isIn;
}

bool Rectangle::isInsideExcludingRectangle(double angleFrom, double angleTo, double x, double y) const {
    // so we have angle range in 0-90, 90-180 or whatever else
    // let's first calculate "excluding rectangle"
    Vector<2> newSizeFrom = Vector<2>();
    Vector<2> newSizeTo = Vector<2>();

    Rectangle rec(*this);
    const double *position = this->getPosition();
    const double *recPosition = rec.getPosition();
    rec.setAngle(angleFrom);
    for (int i = 0; i < 4; i++) {
        newSizeFrom[0] = std::max(newSizeFrom[0], rec.xs[i] - recPosition[0]);
        newSizeFrom[1] = std::max(newSizeFrom[1], rec.ys[i] - recPosition[1]);
    }
    rec.setAngle(angleTo);
    for (int i = 0; i < 4; i++) {
        newSizeTo[0] = std::max(newSizeTo[0], rec.xs[i] - recPosition[0]);
        newSizeTo[1] = std::max(newSizeTo[1], rec.ys[i] - recPosition[1]);
    }
    double newHalfA = std::min(newSizeFrom[0], newSizeTo[0]) + halfA;
    double newHalfB = std::min(newSizeFrom[1], newSizeTo[1]) + halfB;

    // if the point is outside rectangle, all work is done
    // std::abs takes an int as a parameter so, did you mean std::fabs?
    // it's math abs
    if (abs(x - position[0]) <= newHalfA && abs(y - position[1]) <= newHalfB) {
        return 1;
    }
    return 0;
}

int Rectangle::isInsideCircle(const Vector<2, double> &point, double cx, double cy, double r) const {
    if ((cx - point[0]) * (cx - point[0]) + (cy - point[1]) * (cy - point[1]) <= r * r) {
        return 1;
    }
    return 0;
}

int Rectangle::isInsideRect(const Vector<2, double> &point, double cx, double cy, double aLength, double bLength) const {
    if (abs(cx - point[0]) <= aLength / 2 && abs(cy - point[1]) <= bLength / 2) {
        return 1;
    }
    return 0;
}

Vector<2> Rectangle::rotatePoint(const Vector<2> &point, const Matrix<2, 2> &rotation) const {
    return rotation * point;
}

Vector<2> Rectangle::rotatePoint(const Vector<2> &point, const Matrix<2, 2> &rotation, const Vector<2> &center) const {
    Vector<2> translatedPoint = Vector<2>();
    translatedPoint[0] = point[0] - center[0];
    translatedPoint[1] = point[1] - center[1];
    Vector<2> rotated = rotation * translatedPoint;
    rotated[0] += center[0];
    rotated[1] += center[1];
    return rotated;
}

bool Rectangle::isIntersection(double p0_x, double p0_y, double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y) const {
    double s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;
    s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;
    s2_y = p3_y - p2_y;

    double s, t;
    s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
        return 1;
    }

    return 0; // No collision
}

int Rectangle::overlap(BoundaryConditions *bc, Shape<2, 1> *s) const {
    Rectangle *other = (Rectangle *) s;
    Rectangle rectangle(*other);
    this->applyBC(bc, &rectangle);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (isIntersection(xs[i], ys[i], xs[i+1], ys[i+1], rectangle.xs[j], rectangle.ys[j], rectangle.xs[j+1], rectangle.ys[j+1])) {
                return 1;
            }
        }
    }

    return 0;
}

std::string Rectangle::toWolfram() const {
    std::stringstream out;

    Vector<2> thisPos(this->getPosition());
    out << std::fixed;
    out << "GeometricTransformation[Rectangle[{" << -halfA << "," << -halfB << "}, {" << halfA << "," << halfB << "}],"
        << std::endl;
    out << "    {RotationMatrix[" << this->getAngle() << "], " << thisPos << "}]";

    return out.str();
}


std::string Rectangle::toString() const {
    std::stringstream out;

    Vector<2> thisPos(this->getPosition());
    out << "Rectangle{position: " << thisPos << "; a: " << this->a;
    out << "; b: " << this->b << "; angle: " << this->getAngle() << "}";
    return out.str();
}

std::string Rectangle::toPovray() const {
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);
	out << "  polygon{5, <0.5, 0.5, 0.0002>, <-0.5, 0.5, 0.0002>, <-0.5, -0.5, 0.0002>, <0.5, -0.5, 0.0002>, <0.5, 0.5, 0.0002>" << std::endl;;
	out << "	scale <" << this->a << ", " << this->b << ", 1.0>" << std::endl;
	out << "	rotate <0, 0, " << (180*this->getAngle()/M_PI) << ">" << std::endl;
	out << "	translate <";
	for(unsigned short i=0; i<2; i++)
		out << (this->getPosition()[i]) << ", ";
	out << "0.0>" << std::endl << "	texture { pigment { color Red } }" << std::endl << "}" << std::endl;

	return out.str();
}

void Rectangle::store(std::ostream &f) const {
    Shape::store(f);
    f.write((char *) (&this->a), sizeof(double));
    f.write((char *) (&this->b), sizeof(double));
}

void Rectangle::restore(std::istream &f) {
    Shape::restore(f);
    f.read((char *) (&this->a), sizeof(double));
    f.read((char *) (&this->b), sizeof(double));
    this->halfA = a / 2;
    this->halfB = b / 2;
    this->setAngle(this->getAngle());
}