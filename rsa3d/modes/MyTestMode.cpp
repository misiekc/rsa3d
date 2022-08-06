//
// Created by pkua on 06.11.2021.
//

#include <sstream>
#include <memory>
#include <utility>

#include "MyTestMode.h"
#include "../shape/shapes/RoundedCone.h"
#include "../shape/ShapeFactory.h"
#include "../shape/Shape.h"
#include "../RND.h"
#include "../geometry/Vector.h"
#include "../boundary_conditions/FreeBC.h"

void MyTestMode::run() {
	RND rnd;
    FreeBC<3> bc;
	RoundedCone rc1(Matrix<3,3>::identity());
	RoundedCone rc2(Matrix<3,3>::identity());
//	rc1.translate(Vector<3>({2.2075995719042707, 6.0603573519848641, 8.3910673826692808}));
//	rc2.translate(Vector<3>({2.198574832149657, 6.3546092463432187, 8.1357708437132725}));
//	std::cout << (rc1.getPosition()-rc2.getPosition()).norm() << ": " << rc1.overlap(&bc, &rc2) << std::endl;


	rc1.translate(-rc1.getPosition());
	rc2.translate(-rc2.getPosition());
	double dx=0.01;
	for(int i=1; i<=1.0/dx; i++){
		bool b = rc1.overlap(&bc, &rc2);
		std::cout << (rc1.getPosition()-rc2.getPosition()).norm() << ": " << b << std::endl;
		rc2.translate(Vector<3>({0, 0, dx}));
	}

}

void MyTestMode::printHelp(std::ostream &out, const ProgramArguments &arguments) {
    out << arguments.formatUsage("") << std::endl;
    out << std::endl;
    out << "My own test mode for various purposes" << std::endl;
}
