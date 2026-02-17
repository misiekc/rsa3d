/*
* Created by Michal Ciesla on 30.06.2022.
*
* Author: Michal Ciesla & Agnieszka Rejek
*
*/

#include "Superellipsoid.h"
#include "../../utils/Assertions.h"
#include "../OrderParameters.h"

#include <cmath>

double Superellipsoid::a{};
double Superellipsoid::b{};
double Superellipsoid::c{};
double Superellipsoid::p{};

void Superellipsoid::initClass(const std::string &attr) {
    std::istringstream attrStream(attr);

    attrStream >> a >> b >> c >> p;
    if (!attrStream)    throw ValidationException("Wrong attr format");
    ValidateMsg(a > 0 && b > 0 && c > 0, "Superellipsoid semiaxes should be positive");
    ValidateMsg(p > 0, "exponent p should be positive");

    // Axes of increasing size a <= b <= c - the largest axis lies on z-axis
    if (a > b) std::swap(a, b);
    if (a > c) std::swap(a, c);
    if (b > c) std::swap(b, c);

    normalizeVolume();

    ShapeStaticInfo<3, 0> shapeInfo;
    shapeInfo.setCircumsphereRadius(std::cbrt(c*c + b*b + a*a));
    shapeInfo.setInsphereRadius(a*b*c / std::sqrt(a*a*b*b + a*a*c*c + b*b*c*c));
    shapeInfo.setExclusionZoneMaxSpan(0.1*a);
    shapeInfo.setCreateShapeImpl([](RND *rnd) -> Shape* {
        return new Superellipsoid(Matrix<3, 3>::rotation(
                2 * M_PI * rnd->nextValue(),
                std::asin(2 * rnd->nextValue() - 1),
                2 * M_PI * rnd->nextValue())
                );
    });

    Shape::setShapeStaticInfo(shapeInfo);
}

void Superellipsoid::normalizeVolume() {
    double volume = (2.0/3.0)*a*b*c/(p*p);

//    volume *= std::beta(1.0/(2.0*p), 1.0/(2.0*p));
//    volume *= std::beta(1.0/p, 1.0/(2.0*p));
	double factor = 1/std::cbrt(volume);
	double gamma_1_p = std::tgamma(1.0/p);
	double gamma_1_2p = std::tgamma(1.0/(2.0*p));
	volume *= gamma_1_2p*gamma_1_2p / gamma_1_p;
	volume *= gamma_1_p*gamma_1_2p / std::tgamma(3.0/(2.0*p));

    Superellipsoid::a *= factor;
    Superellipsoid::b *= factor;
    Superellipsoid::c *= factor;
}

Superellipsoid::Superellipsoid(const Matrix<3, 3> &orientation)
{
	this->orientation = orientation;
	this->orientationTr = this->orientation.transpose();
	this->orientationInv = this->orientation.inverse();

}

// r0 position

Vector<3> Superellipsoid::r_relative(const Vector<3> &r) const{
	return (this->orientationInv)*(r - this->getPosition());
}

double Superellipsoid::inner_shape_function(const Vector<3> &r) const{
	Vector<3> r_relative = this->r_relative(r);
	return 	std::pow(std::abs(r_relative[0] / a), 2*p) +
			std::pow(std::abs(r_relative[1] / b), 2*p) +
			std::pow(std::abs(r_relative[2] / c), 2*p);
}

double Superellipsoid::shape_function(const Vector<3> &r) const{
	return std::pow(this->inner_shape_function(r), 1.0/Superellipsoid::p) - 1;
}

Vector<3> Superellipsoid::inner_shape_function_gradient(const Vector<3> &r) const{
	Vector<3> gradient;
	Vector<3> r_relative = this->r_relative(r);
	gradient[0] = (2*p/a)*std::pow(std::abs(r_relative[0] / a), 2*p - 1);
	gradient[1] = (2*p/b)*std::pow(std::abs(r_relative[1] / b), 2*p - 1);
	gradient[2] = (2*p/c)*std::pow(std::abs(r_relative[2] / c), 2*p - 1);
	return this->orientation*gradient;
}

Vector<3> Superellipsoid::shape_function_gradient(const Vector<3> &r) const{
	Vector<3> gradient;
	double factor = (1.0/p)*std::pow(this->inner_shape_function(r), 1.0/p - 1);
	gradient = factor*this->inner_shape_function_gradient(r);
	return gradient;
}

Matrix<3, 3> Superellipsoid::inner_shape_function_hesian(const Vector<3> &r) const{
	Matrix<3, 3> hesian{};
	Vector<3> r_relative = this->r_relative(r);
	hesian(0, 0) = (2*p*(2*p-1)*std::pow(std::abs(r_relative[0]/a), 2*p-2))/(a*a);
	hesian(0, 1) = 0.0;
	hesian(0, 2) = 0.0;
	hesian(1, 0) = 0.0;
	hesian(1, 1) = (2*p*(2*p-1)*std::pow(std::abs(r_relative[1]/b), 2*p-2))/(b*b);
	hesian(1, 2) = 0.0;
	hesian(2, 0) = 0.0;
	hesian(2, 1) = 0.0;
	hesian(2, 2) = (2*p*(2*p-1)*std::pow(std::abs(r_relative[2]/c), 2*p-2))/(c*c);
	return this->orientation*hesian*this->orientationTr;
}

Matrix<3, 3> Superellipsoid::shape_function_hesian(const Vector<3> &r) const{
	Matrix<3, 3> hesian;
	double inner_shape = this->inner_shape_function(r);
	Vector<3> inner_gradient = this->inner_shape_function_gradient(r);
	Matrix<3, 3> inner_hesian = this->inner_shape_function_hesian(r);

	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			hesian(i,j) = 	(1/p)*
								((1.0/p - 1)*std::pow(inner_shape, 1.0/p - 2)*inner_gradient[i]*inner_gradient[j] +
								std::pow(inner_shape, 1.0/p - 1)*inner_hesian(i,j));
	return hesian;
}

bool Superellipsoid::overlap(BoundaryConditions<3> *bc, const Shape<3, 0> *s) const {
    	switch (this->overlapEarlyRejection(bc, s)) {
    		case TRUE:      return true;
    		case FALSE:     return false;
    		case UNKNOWN:   break;
    	}
    	Superellipsoid superellipsoid = dynamic_cast<const Superellipsoid &>(*s);
    	this->applyBC(bc, &superellipsoid);


    	double lambda = 0.5;
    	Vector<3> rc = (this->getPosition() + superellipsoid.getPosition())*0.5;
    	double delta_lambda = 0;
    	Vector<3> delta_rc;
    	delta_rc[0] = 0.0; delta_rc[1] = 0.0; delta_rc[2] = 0.0;

    	double pw_potential;
    	for(int i=0; i<20; i++){
    		lambda = lambda + delta_lambda;
    		rc = rc + delta_rc;
/*
    		std::cout<< "--------" << i << "--------" << std::endl;
    		std::cout<< "lambda: " << lambda << std::endl;
    		std::cout<< "delta_lambda: " << delta_lambda << std::endl;
    		std::cout<< "rc: " << rc << std::endl;
    		std::cout<< "delta_rc: " << delta_rc << std::endl;
    		std::cout<< "shape_functions (rc): (" << this->shape_function(rc) << ", " << superellipsoid.shape_function(rc) << ")" << std::endl;
*/

    		Vector<3> nabla_zeta_A = this->shape_function_gradient(rc);
//    		std::cout<< "nabla_A: " << nabla_zeta_A << std::endl;
    		Vector<3> nabla_zeta_B = superellipsoid.shape_function_gradient(rc);
//    		std::cout<< "nabla_B: " << nabla_zeta_B << std::endl;
    		Vector<3> delta_g = nabla_zeta_B - nabla_zeta_A;
//       	std::cout<< "delta_g: " << delta_g << std::endl;

    		Matrix<3, 3> matrix_M = this->shape_function_hesian(rc)*lambda + superellipsoid.shape_function_hesian(rc)*(1 - lambda);
//       	std::cout<< "M: " << matrix_M << std::endl;
    		Matrix<3, 3> matrix_MInv = matrix_M.inverse();
    		double zeta_ll = delta_g*(matrix_MInv*delta_g);

    		Vector<3> nabla_zeta_AB = nabla_zeta_A*lambda + nabla_zeta_B*(1-lambda);

    		delta_lambda = (-1.0/zeta_ll)*((superellipsoid.shape_function(rc) - this->shape_function(rc))  - delta_g*(matrix_MInv*nabla_zeta_AB));

    		delta_rc = matrix_MInv*(delta_g*delta_lambda - nabla_zeta_AB);

    		pw_potential = this->shape_function(rc)*lambda + superellipsoid.shape_function(rc)*(1-lambda);
//    		std::cout << "potential: " << pw_potential << std::endl;
     	}
    	return pw_potential<0;
}

bool Superellipsoid::pointInside(BoundaryConditions<3> *bc, const Vector<3> &position, const Orientation<0> &orientation, const Orientation<0> &orientationRange) const {
    // Transform point coordinates to Ellipsoid coordinate system
    Vector<3> bcPos = position + bc->getTranslation(this->getPosition(), position);
    Vector<3> pointAligned = this->orientationTr * (bcPos - this->getPosition());

    double d = std::pow(pointAligned[0]/Superellipsoid::a, 2*Superellipsoid::p) + std::pow(pointAligned[1]/Superellipsoid::b, 2*Superellipsoid::p) + std::pow(pointAligned[2]/Superellipsoid::c, 2*Superellipsoid::p);
    return (d <= 1.0 + this->getShapeStaticInfo().getInsphereRadius());
}

void Superellipsoid::store(std::ostream &f) const {
    Shape::store(f);
    double d;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
            d = this->orientation(i, j);
            f.write((char *)(&d), sizeof(double));
        }
    }
}

void Superellipsoid::restore(std::istream &f) {
    Shape::restore(f);
    double d;
    for (unsigned short i=0; i<3; i++){
        for (unsigned short j=0; j<3; j++){
            f.read((char *)&d, sizeof(double));
            this->orientation(i, j) = d;
        }
    }
	this->orientationTr = this->orientation.transpose();
	if (this->orientation.det()!=0)
		this->orientationInv = this->orientation.inverse();
}

Shape<3, 0> *Superellipsoid::clone() const {
    return new Superellipsoid(*this);
}

std::string Superellipsoid::toPovray() const {
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);
    Vector<3> position = this->getPosition();
	out << "  superellipsoid { <" << 1.0/p << ", " << 1.0/p << ">" << std::endl;
	out << "    scale < " << a << ", " << b << ", " << c << " >" << std::endl;
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

std::string Superellipsoid::toWolfram() const {
    std::stringstream out;
    out << std::fixed;
    out << "GeometricTransformation[" << std::endl;
    out << "    Ellipsoid[{0, 0, 0}, {" << a << ", " << b << ", " << c << "}]," << std::endl;
    out << "    AffineTransform[" << std::endl;
    out << "        {{{" << this->orientation(0, 0) << ", " << this->orientation(0, 1) << ", " << this->orientation(0, 2) << "}," << std::endl;
    out << "          {" << this->orientation(1, 0) << ", " << this->orientation(1, 1) << ", " << this->orientation(1, 2) << "}," << std::endl;
    out << "          {" << this->orientation(2, 0) << ", " << this->orientation(2, 1) << ", " << this->orientation(2, 2) << "}}," << std::endl;
    out << "          " << this->getPosition() << "}]]";
    return out.str();
}

std::vector<double> Superellipsoid::calculateOrder(const OrderCalculable *other) const {
    auto &otherEllipsoid = dynamic_cast<const Superellipsoid&>(*other);    // use reference so dynamic_cast would throw on error
    Matrix<3,3> product = this->orientation.transpose() * otherEllipsoid.orientation;

    double v = product(0,0);
    if(a==b){
    	v = product(2,2);
    }

    std::vector<double> result(5);
    result[0] = integer_legendre::P2(v);
    result[1] = (product(0, 0) + product(1, 1) + product(2, 2));
    return result;
}



