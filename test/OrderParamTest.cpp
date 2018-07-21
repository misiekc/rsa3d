//
// Created by PKua on 21.07.18.
//

#include <cstdlib>
#include <memory>

#include "OrderParamTest.h"
#include "../rsa3d/Utils.h"
#include "../rsa3d/shape/ShapeFactory.h"
#include "../rsa3d/shape/shapes/regular_solid/RegularSolidBase.h"
#include "../rsa3d/shape/OrderCalculable.h"

namespace {
    using shape_ptr = std::unique_ptr<RSAShape>;
    using shape_cptr = std::unique_ptr<const RSAShape>;

    struct Boundaries {
        double      minParamValue = std::numeric_limits<double>::infinity();
        double      maxParamValue = -std::numeric_limits<double>::infinity();
        shape_cptr  minShape;
        shape_cptr  maxShape;
    };

    shape_ptr acquire_first_shape(const std::string &particleName) {
        ShapeFactory::initShapeClass(particleName, "");
        RND rnd;
        return shape_ptr(ShapeFactory::createShape(&rnd));
    }

    bool is_shape_compatible(const std::string &shapeName, const RSAShape *shape) {
        auto regularSolid = dynamic_cast<const RegularSolidBase*>(shape);
        auto orderCalculable = dynamic_cast<const OrderCalculable*>(shape);
        return regularSolid != nullptr && orderCalculable != nullptr;
    }

    double calculate_full_order(const RSAShape &firstShape, const RSAShape &secondShape) {
        auto &firstOrderCalc = dynamic_cast<const OrderCalculable&>(firstShape);
        auto &secondOrderCalc = dynamic_cast<const OrderCalculable&>(secondShape);
        auto orders = firstOrderCalc.calculateOrder(&secondOrderCalc);

        return orders[1];
    }

    Boundaries find_boundaries(const RSAShape &firstShape, std::size_t iterations) {
        RND rnd{};
        Boundaries result{};

        std::cout << "[INFO] Finding boundaries... " << std::flush;
        for (std::size_t i = 0; i < iterations; i++) {
            auto secondShape = shape_ptr(ShapeFactory::createShape(&rnd));

            double fullOrder = calculate_full_order(firstShape, *secondShape);
            if (fullOrder < result.minParamValue) {
                result.minParamValue = fullOrder;
                result.minShape.reset(secondShape->clone());
            }
            if (fullOrder > result.maxParamValue) {
                result.maxParamValue = fullOrder;
                result.maxShape.reset(secondShape->clone());
            }
        }

        std::cout << "DONE" << std::endl;
        return result;
    }

    double find_disorder_value(const RSAShape &firstShape, std::size_t iterations) {
        RND rnd;
        double result = 0;

        std::cout << "[INFO] Finding disorder value... " << std::flush;
        for (std::size_t i = 0; i < iterations; i++) {
            auto secondShape = shape_ptr(ShapeFactory::createShape(&rnd));
            double fullOrder = calculate_full_order(firstShape, *secondShape);
            result += fullOrder;
        }

        std::cout << "DONE" << std::endl;
        return result / iterations;
    }

    Matrix<3, 3> rotation_around_axis(const Vector<3> &axis, double angle) {
        double sina = std::sin(angle);
        double cosa = std::cos(angle);
        Matrix<3, 3> K = {{      0 , -axis[2],  axis[1],
                            axis[2],       0,  -axis[0],
                           -axis[1],  axis[0],       0 }};

        // Rodrigues' rotation formula
        return Matrix<3, 3>::identity() + K*sina + (1 - cosa)*K*K;
    }

    double find_nematic_order_value(const RSAShape &firstShape, std::size_t iterations) {
        RND rnd;
        auto &firstRegular = dynamic_cast<const RegularSolidBase&>(firstShape);
        Vector<3> fixedAxis = firstRegular.getFaceAxes()[0];

        double result = 0;
        std::cout << "[INFO] Finding nematic order value... " << std::flush;
        for (std::size_t i = 0; i < iterations; i++) {
            auto secondShape = shape_ptr(firstShape.clone());
            auto &secondRegular = dynamic_cast<RegularSolidBase&>(*secondShape);

            Matrix<3, 3> oldOrient = secondRegular.getOrientationMatrix();
            double randomAngle = 2*M_PI*rnd.nextValue();
            Matrix<3, 3> newOrient = rotation_around_axis(fixedAxis, randomAngle) * oldOrient;
            secondRegular.setOrientationMatrix(newOrient);

            double fullOrder = calculate_full_order(firstShape, *secondShape);
            result += fullOrder;
        }

        std::cout << "DONE" << std::endl;
        return result / iterations;
    }

    void report_results(const RSAShape &first, const Boundaries &boundaries, double disorder, double nematicOrder) {
        std::cout << std::endl;
        std::cout << "[INFO] Min value           : " << boundaries.minParamValue << std::endl;
        std::cout << "[INFO] Max value           : " << boundaries.maxParamValue << std::endl;
        std::cout << "[INFO] Disorder value      : " << disorder << std::endl;
        std::cout << "[INFO] Nematic order value : " << nematicOrder << std::endl;
    }
}

int order_param_test::main(int argc, char **argv) {
    if (argc < 4)
        die("[ERROR] Usage: ./rsa_test order_param_test [particle] [iterations]");

    auto shape = acquire_first_shape(argv[2]);
    if (!is_shape_compatible(argv[2], shape.get()))
        die("[ERROR] " + std::string(argv[2]) + " is not both RegularSolid and OrderCalculable");

    std::size_t iterations = std::stoul(argv[3]);
    if (iterations == 0)
        die("[ERROR] iterations == 0");

    Boundaries boundaries = find_boundaries(*shape, iterations);
    double disorder = find_disorder_value(*shape, iterations);
    double nematicOrder = find_nematic_order_value(*shape, iterations);

    report_results(*shape, boundaries, disorder, nematicOrder);
    return EXIT_SUCCESS;
}