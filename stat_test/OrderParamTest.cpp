//
// Created by PKua on 21.07.18.
//

#include <cstdlib>
#include <memory>
#include <tuple>

#include "OrderParamTest.h"
#include "../rsa3d/shape/ShapeFactory.h"
#include "../rsa3d/shape/shapes/regular_solid/RegularSolidBase.h"
#include "utils/ShapeGenerators.h"
#include "utils/ShapePairFactory.h"
#include "../rsa3d/boundary_conditions/FreeBC.h"
#include "utils/InfoLooper.h"
#include "../rsa3d/utils/Quantity.h"
#include "utils/TestExitCodes.h"

namespace {
    using shape_ptr = std::unique_ptr<RSAShape>;
    using shape_cptr = std::unique_ptr<const RSAShape>;

    struct Boundaries {
        double      minParamValue = std::numeric_limits<double>::infinity();
        double      maxParamValue = -std::numeric_limits<double>::infinity();
        shape_cptr  minShape;
        shape_cptr  maxShape;
    };

    struct Result {
        Boundaries boundaries;
        Quantity disorder;
        Quantity nematicOrder;
        Quantity nematicOrderFullRotations;
        int nematicExitCode;
        std::unique_ptr<RSAShape> firstShape;

        void report(std::ostream &out) const;

        int getExitCode();
    };

    bool is_shape_compatible(const RSAShape *shape) {
        auto regularSolid = dynamic_cast<const RegularSolidBase*>(shape);
        auto orderCalculable = dynamic_cast<const OrderCalculable*>(shape);
        return regularSolid != nullptr && orderCalculable != nullptr;
    }

    double calculate_full_order(const RSAShape &firstShape, const RSAShape &secondShape) {
        auto &firstOrderCalc = dynamic_cast<const OrderCalculable&>(firstShape);
        auto &secondOrderCalc = dynamic_cast<const OrderCalculable&>(secondShape);
        auto orders = firstOrderCalc.calculateOrder(&secondOrderCalc);

        return orders[0];
    }

    Boundaries find_boundaries(const RSAShape &firstShape, std::size_t iterations) {
        RND rnd{};
        Boundaries result{};

        std::cout << "[INFO] Finding boundaries... " << std::flush;
        for (std::size_t i = 0; i < iterations; i++) {
            auto secondShape = generate_randomly_oriented_shape(Vector<3>(), &rnd);

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

    Quantity find_disorder_value(const RSAShape &firstShape, std::size_t iterations) {
        RND rnd;

        std::vector<double> values;
        std::cout << "[INFO] Finding disorder value... " << std::flush;
        for (std::size_t i = 0; i < iterations; i++) {
            auto secondShape = generate_randomly_oriented_shape(Vector<3>(), &rnd);
            double fullOrder = calculate_full_order(firstShape, *secondShape);
            values.push_back(fullOrder);
        }

        std::cout << "DONE" << std::endl;

        Quantity result;
        result.separator = Quantity::PLUS_MINUS;
        result.calculateFromSamples(values);
        return result;
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

    Quantity find_nematic_order_value(const RSAShape &firstShape, std::size_t iterations) {
        RND rnd;
        auto &firstRegular = dynamic_cast<const RegularSolidBase&>(firstShape);
        Vector<3> fixedAxis = firstRegular.getFaceAxes()[0];

        std::vector<double> values;
        std::cout << "[INFO] Finding nematic order value... " << std::flush;
        for (std::size_t i = 0; i < iterations; i++) {
            auto secondShape = shape_ptr(firstShape.clone());
            auto &secondRegular = dynamic_cast<RegularSolidBase&>(*secondShape);

            Matrix<3, 3> oldOrient = secondRegular.getOrientationMatrix();
            double randomAngle = 2*M_PI*rnd.nextValue();
            Matrix<3, 3> newOrient = rotation_around_axis(fixedAxis, randomAngle) * oldOrient;
            secondRegular.setOrientationMatrix(newOrient);

            double fullOrder = calculate_full_order(firstShape, *secondShape);
            values.push_back(fullOrder);
        }

        std::cout << "DONE" << std::endl;

        Quantity result;
        result.separator = Quantity::PLUS_MINUS;
        result.calculateFromSamples(values);
        return result;
    }

    /* Averages order parameter over firstShape and a bunch of second shapes translated in the direction of fixedAxis
     * just a bit more than the insphere radius, however there much be no overlap. */
    auto find_nematic_order_value_for_axis(const RSAShape &firstShape, size_t iterations, const Vector<3> &fixedAxis) {
        // Shapes will be separated by a little bit more than the minimal distance
        Vector<3> translation = fixedAxis * (1.01 * 2 * RSAShape::getInsphereRadius());

        double numberOfNonoverlaping{};
        RND rnd;
        RSAFreeBC bc;
        InfoLooper looper(iterations * 100, iterations, "pairs tested...");
        std::vector<double> values;
        while(looper.step()) {
            auto secondShape = generate_randomly_oriented_shape(Vector<3>(), &rnd);
            secondShape->translate(translation);
            if (firstShape.overlap(&bc, secondShape.get()))
                continue;

            numberOfNonoverlaping++;
            double fullOrder = calculate_full_order(firstShape, *secondShape);
            values.push_back(fullOrder);
        }

        Quantity result;
        result.separator = Quantity::PLUS_MINUS;
        result.calculateFromSamples(values);
        return std::make_tuple(numberOfNonoverlaping, result);
    }

    /* Returns a tuple: test exit code and the average value of order parameter for solids almost closest possible
     * Test exit code is TEST_ERROR if all sampled solids overlapped - it check face axes and vertex axes of symmetry
     * Platonic solid */
    auto find_nematic_order_value_full_rotations(const RSAShape &firstShape, std::size_t iterations) {
        auto &firstRegular = dynamic_cast<const RegularSolidBase&>(firstShape);
        auto firstSymmetry = firstRegular.createSymmetryShape();

        std::cout << "[INFO] Finding nematic order value for face axis using full rotations... " << std::endl;
        auto [numberOfIterations, result] = find_nematic_order_value_for_axis(firstShape, iterations,
                                                                              firstSymmetry->getFaceAxes()[0]);

        if (numberOfIterations != 0) {
            std::cout << "DONE. Iterations without overlap: " << numberOfIterations << std::endl;
            return std::make_tuple(TEST_SUCCESS, result);
        }

        std::cout << "[INFO] No non-overlapping particles found. Trying vertex axes..." << std::endl;
        std::tie(numberOfIterations, result) = find_nematic_order_value_for_axis(firstShape, iterations,
                                                                                 firstSymmetry->getVertexAxes()[0]);

        if (numberOfIterations != 0) {
            std::cout << "[INFO] DONE. Iterations without overlap: " << numberOfIterations << std::endl;
            return std::make_tuple(TEST_SUCCESS, result);
        } else {
            std::cout << "[ERROR] No non-overlapping particles found. " << std::endl;
            return std::make_tuple(TEST_ERROR, result);
        }
    }

    void Result::report(std::ostream &out) const {
        out << std::endl;
        out << "[INFO] Min value                                : " << this->boundaries.minParamValue << std::endl;
        out << "[INFO] Max value                                : " << this->boundaries.maxParamValue << std::endl;
        out << "[INFO] Disorder value                           : " << this->disorder << std::endl;
        out << "[INFO] Nematic order value                      : " << this->nematicOrder << std::endl;
        if (this->nematicExitCode == TEST_SUCCESS)
            out << "[INFO] Nematic order value using full rotations : " << this->nematicOrderFullRotations << std::endl;
        else
            out << "[INFO] Nematic order value using full rotations : ERROR" << std::endl;
        out << std::endl;
        out << "[INFO] Configuration yielding minimal value:" << std::endl;

        RSAShapePairFactory::ShapePair pair = {this->firstShape->clone(), this->boundaries.minShape->clone()};
        pair.print(std::cout);
    }

    int Result::getExitCode() {
        return this->nematicExitCode;
    }

    Result perform_tests(std::unique_ptr<RSAShape> firstShape, size_t iterations) {
        Result result;
        result.firstShape = std::move(firstShape);
        result.boundaries = find_boundaries(*result.firstShape, iterations);
        result.disorder = find_disorder_value(*result.firstShape, iterations);
        result.nematicOrder = find_nematic_order_value(*result.firstShape, iterations);
        std::tie(result.nematicExitCode, result.nematicOrderFullRotations)
                = find_nematic_order_value_full_rotations(*result.firstShape, iterations);
        return result;
    }
}

int order_param_test::main(int argc, char **argv) {
    if (argc < 5) {
        std::cerr << "[ERROR] Usage: ./rsa_test order_param_test [particle] [attributes] [iterations]" << std::endl;
        return TEST_ERROR;
    }

    std::string particleName = argv[2];
    std::string particleAttr = argv[3];
    ShapeFactory::initShapeClass(particleName, particleAttr);

    RND rnd;
    auto firstShape = generate_randomly_oriented_shape(Vector<3>(), &rnd);
    if (!is_shape_compatible(firstShape.get())) {
        std::cerr << "[ERROR] " << particleName << " is not both RegularSolid and OrderCalculable" << std::endl;
        return TEST_ERROR;
    }

    std::size_t iterations = std::stoul(argv[4]);
    if (iterations == 0) {
        std::cerr << "[ERROR] iterations == 0" << std::endl;
        return TEST_ERROR;
    }

    Result result = perform_tests(std::move(firstShape), iterations);
    result.report(std::cout);
    return result.getExitCode();
}
