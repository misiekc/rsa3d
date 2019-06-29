/*
* Created by PKua on 17.07.18.
*
* Author: Wojciech Nowak
*
*/

#include "Ellipsoid.h"
#include "../../utils/Assertions.h"
#include "../OrderParameters.h"

double Ellipsoid::a{};
double Ellipsoid::b{};
double Ellipsoid::c{};

void Ellipsoid::initClass(const std::string &attr) {
    std::istringstream attrStream(attr);

    attrStream >> a >> b >> c;
    if (!attrStream)    throw ValidationException("Wrong attr format");
    ValidateMsg(a > 0 && b > 0 && c > 0, "Ellipsoid semiaxes should be positive");

    // Axes of increasing size a <= b <= c - the largest axis lies on z-axis
    if (a > b) std::swap(a, b);
    if (a > c) std::swap(a, c);
    if (b > c) std::swap(b, c);

    normalizeVolume();

    ShapeStaticInfo<3, 0> shapeInfo;
    shapeInfo.setCircumsphereRadius(c);
    shapeInfo.setInsphereRadius(a);
    shapeInfo.setCreateShapeImpl([](RND *rnd) -> Shape* {
        return new Ellipsoid(Matrix<3, 3>::rotation(
                2 * M_PI * rnd->nextValue(),
                std::asin(2 * rnd->nextValue() - 1),
                2 * M_PI * rnd->nextValue()));
    });

    Shape::setShapeStaticInfo(shapeInfo);
}

void Ellipsoid::normalizeVolume() {
    double volume = 4./3 * M_PI * a * b * c;
    double factor = 1/std::cbrt(volume);
    a *= factor;
    b *= factor;
    c *= factor;
}

Ellipsoid::Ellipsoid(const Matrix<3, 3> &orientation) : orientation(orientation), X{{a*a,0,0,0,b*b,0,0,0,c*c}},
                                                        M(matrixM())
{

}

bool Ellipsoid::overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const {
    Ellipsoid ellipsoid = dynamic_cast<const Ellipsoid &>(*s);
    this->applyBC(bc, &ellipsoid);

    Vector<3> rAB = vectorRAB(this->getPosition(), ellipsoid.getPosition());

    // Early rejection
    double rABNorm2 = rAB.norm2();
    if (rABNorm2 > 4*c*c)
        return false;
    else if (rABNorm2 < 4*a*a)
        return true;

    double lambda = INITIAL_GUESS;
	double fder = 0.0;
	double sder = 0.0;
	double derivativeQuotient = 0.0;
	double fAB = 0.0;
	double previousLambda = 0.0;

	for (int i = 0; i < 20; i++) {
		try {
			fder = firstDerivative(this->getEllipsoidMatrix(), ellipsoid.getEllipsoidMatrix(), rAB, lambda);
            sder = secondDerivative(this->getEllipsoidMatrix(), ellipsoid.getEllipsoidMatrix(), rAB, lambda);

			derivativeQuotient = fder / sder;
			previousLambda = lambda;
			lambda = getNewLambda(lambda, derivativeQuotient);
		}
		catch (std::overflow_error e) {
			std::cout << e.what();
		}
		if (std::abs(fder) < 1e-8 || std::abs(lambda - previousLambda) < 1e-8 ) {
            fAB = overlapFunction(this->getEllipsoidMatrix(), ellipsoid.getEllipsoidMatrix(), rAB, lambda);

			if (fAB < 1) {
				return true;    // overlap
			}
			else {
				return false;   // disjoint
			}
		}
	}
	return false;
}

bool Ellipsoid::pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation, double orientationRange) const {
    // Transform point coordinates to Ellipsoid coordinate system
    Vector<3> bcPos = position + bc->getTranslation(this->getPosition(), position);
    Vector<3> pointAligned = this->orientation.transpose() * (bcPos - this->getPosition());

    double d = pointAligned[0]*pointAligned[0]/(4*Ellipsoid::a*Ellipsoid::a) + pointAligned[1]*pointAligned[1]/((Ellipsoid::b + Ellipsoid::a)*(Ellipsoid::b + Ellipsoid::a)) + pointAligned[2]*pointAligned[2]/((Ellipsoid::c + Ellipsoid::a)*(Ellipsoid::c + Ellipsoid::a));
    return (d <= 1.0);
}

void Ellipsoid::store(std::ostream &f) const {
    Shape::store(f);
    double d;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
            d = this->orientation(i, j);
            f.write((char *)(&d), sizeof(double));
        }
    }
}

void Ellipsoid::restore(std::istream &f) {
    Shape::restore(f);
    double d;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
            f.read((char *)&d, sizeof(double));
            this->orientation(i, j) = d;
        }
    }
    this->M = this->matrixM();
}

Shape<3, 0> *Ellipsoid::clone() const {
    return new Ellipsoid(*this);
}

std::string Ellipsoid::toPovray() const {
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);
    Vector<3> position = this->getPosition();
	out << "  sphere { < 0, 0, 0 >, 1.0" << std::endl;
	out << "    scale < " << Ellipsoid::a << ", " << Ellipsoid::b << ", " << Ellipsoid::c << " >" << std::endl;
	out << "    matrix <" << std::endl;
	out << "      " << this->orientation(0, 0) << ", " << this->orientation(1, 0) << ", " << this->orientation(2, 0) << ", " << std::endl;
    out << "      " << this->orientation(0, 1) << ", " << this->orientation(1, 1) << ", " << this->orientation(2, 1) << ", " << std::endl;
    out << "      " << this->orientation(0, 2) << ", " << this->orientation(1, 2) << ", " << this->orientation(2, 2) << ", " << std::endl;
    out << "      " << position[0] << ", " << position[1] << ", " << position[2] << std::endl;
    out << "    >" << std::endl;
	out << "    texture { pigment { color Red } }" << std::endl;
	out << "  }" << std::endl;
	return out.str();
}

std::string Ellipsoid::toWolfram() const {
    std::stringstream out;
    out << std::fixed;
    out << "GeometricTransformation[" << std::endl;
    out << "    Ellipsoid[{0, 0, 0}, {" << Ellipsoid::a << ", " << Ellipsoid::b << ", " << Ellipsoid::c << "}]," << std::endl;
    out << "    AffineTransform[" << std::endl;
    out << "        {{{" << this->orientation(0, 0) << ", " << this->orientation(0, 1) << ", " << this->orientation(0, 2) << "}," << std::endl;
    out << "          {" << this->orientation(1, 0) << ", " << this->orientation(1, 1) << ", " << this->orientation(1, 2) << "}," << std::endl;
    out << "          {" << this->orientation(2, 0) << ", " << this->orientation(2, 1) << ", " << this->orientation(2, 2) << "}}," << std::endl;
    out << "          " << this->getPosition() << "}]]";
    return out.str();
}

Matrix<3, 3> Ellipsoid::getEllipsoidMatrix() const {
    return M;
}


inline Matrix<3, 3> Ellipsoid::matrixM() const {
    return (orientation * X) * orientation.transpose();
}

inline Matrix<3, 3> Ellipsoid::matrixC(const Matrix<3, 3> &A, const Matrix<3, 3> &B,  double lambda) const {
    return ((B * lambda) + (A * (1 - lambda))).inverse();
}

inline double Ellipsoid::quadraticForm(const Matrix<3, 3> &M, const Vector<3> &v) const {
	Vector<3> vM = M * v;
    return v * vM;
}

inline Vector<3> Ellipsoid::vectorRAB(const Vector<3>& rA, const Vector<3>& rB) const {
    return rB - rA;
}

inline double Ellipsoid::overlapFunction(const Matrix<3, 3> &A, const Matrix<3, 3> &B, const Vector<3> &rAB, double lambda) const {
    return lambda * (1 - lambda) * quadraticForm(matrixC(A, B, lambda), rAB);
}

inline double Ellipsoid::firstDerivative(const Matrix<3, 3> &A, const Matrix<3, 3> &B, const Vector<3> &rAB, double lambda) const {
    Matrix<3, 3> M = B - A;
    return ((((1 - 2 * lambda) * ((A + M * lambda).inverse())) + ((lambda*lambda - lambda)*((A + M * lambda).inverse()) * M * ((A + M * lambda).inverse()))).transpose()*rAB)*rAB;
}

inline double Ellipsoid::secondDerivative(const Matrix<3, 3> &A, const Matrix<3, 3> &B, const Vector<3> &rAB, double lambda) const {
	Matrix<3, 3> C = matrixC(A, B, lambda);
	Matrix<3, 3> AC = A * C;
	Matrix<3, 3> CAC = C * AC;
	Matrix<3, 3> BC = B * C;
	Matrix<3, 3> CBC = C * BC;
	Matrix<3, 3> CACBC = CAC * BC;
	Matrix<3, 3> CBCAC = CBC * AC;

	return (-1.0)*(quadraticForm(CACBC, rAB) + quadraticForm(CBCAC, rAB));
 }

double Ellipsoid::getNewLambda(double lambda, double derivativeQuotient) const {
    double temp = lambda - derivativeQuotient;
	if (temp <= 0)
		return lambda / 2.;
	else if (temp >= 1)
		return (lambda + 1.) / 2.;
	else
		return temp;
}

std::vector<double> Ellipsoid::calculateOrder(const OrderCalculable *other) const {
    auto &otherEllipsoid = dynamic_cast<const Ellipsoid&>(*other);    // use reference so dynamic_cast would throw on error
    Matrix<3,3> product = this->orientation.transpose() * otherEllipsoid.orientation;

    double v = product(0,0);
    if(a==b){
    	v = product(2,2);
    }

    std::vector<double> result(5);
    result[0] = legendre::P2(v);
    result[1] = (product(0, 0) + product(1, 1) + product(2, 2));
    return result;
}
