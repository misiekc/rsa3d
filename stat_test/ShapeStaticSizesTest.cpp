//--------------------------------------------------------------------------------------------
// Test of cuboid intersection algorithm. It uses standard triangle-triangle intersection
// algorithm to check results given by Cuboid::overlap
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------


#include <iostream>
#include <cstring>
#include <memory>
#include <utility>

#include "ShapeStaticSizesTest.h"
#include "../rsa3d/shape/ShapeFactory.h"
#include "../rsa3d/utils/Utils.h"
#include "../rsa3d/FreeBC.h"
#include "utils/InfoLooper.h"
#include "utils/ShapeGenerators.h"
#include "utils/TestExitCodes.h"

namespace
{
    using ShapePair = RSAShapePairFactory::ShapePair;

    class Result {
    private:
        RSAFreeBC bc;

        std::size_t voxelConflicts{};
        std::size_t biggerVoxelConflicts{};
        ShapePair voxelConflictExample{};

        std::size_t neighbourListConflicts{};
        std::size_t smallerNeighbourListConflicts{};
        ShapePair neighbourListConflictExample{};

    public:
        void evaluateForVoxel(const ShapePair &pair);
        void evaluateForNeighbourList(const ShapePair &pair);
        void evaluateForBiggerVoxel(const ShapePair &pair);
        void evaluateForSmallerNeighbourList(const ShapePair &pair);

        void print(std::ostream &out) const;
        bool success() const;
    };

    void Result::evaluateForVoxel(const ShapePair &pair) {
        if (!pair.first()->overlap(&bc, pair.second())) {
            if (this->voxelConflicts == 0)
                this->voxelConflictExample = pair;
            this->voxelConflicts++;
        }
    }

    void Result::evaluateForBiggerVoxel(const ShapePair &pair) {
        if (!pair.first()->overlap(&bc, pair.second()))
            this->biggerVoxelConflicts++;
    }

    void Result::evaluateForNeighbourList(const ShapePair &pair) {
        if (pair.first()->overlap(&bc, pair.second())) {
            if (this->neighbourListConflicts == 0)
                this->neighbourListConflictExample = pair;
            this->neighbourListConflicts++;
        }
    }

    void Result::evaluateForSmallerNeighbourList(const ShapePair &pair) {
        if (pair.first()->overlap(&bc, pair.second()))
            this->smallerNeighbourListConflicts++;
    }

    /* Presents result on the out ostream */
    void Result::print(std::ostream &out) const {
        out << std::endl;
        if (voxelConflicts != 0) {
            out << "[FAILED] Voxel size too big, or Shape::overlap incorrect. " << std::endl;
            out << "[INFO] Example configuration of two non-overlapping (according to the overlap() method) shapes:";
            out << std::endl;
            voxelConflictExample.print(out);
            out << std::endl << std::endl;
        } else if (biggerVoxelConflicts == 0) {
            out << "[FAILED] Voxel size too small or Shape::overlap incorrect. TRY UP TO 100000 PAIRS" << std::endl;
        } else {
            out << "[PASSED] Voxel size correct and optimal" << std::endl;
        }

        if (neighbourListConflicts != 0) {
            out << "[FAILED] NeighbourGrid cell to small, or Shape::overlap incorrect" << std::endl;
            out << "[INFO] Example configuration of two overlapping (according to the overlap() method) shapes:";
            out << std::endl;
            neighbourListConflictExample.print(out);
            out << std::endl << std::endl;
        } else if (smallerNeighbourListConflicts == 0) {
            out << "[WARNING] NeighbourGrid cell too big or Shape::overlap incorrect. TRY UP TO 100000 PAIRS" << std::endl;
        } else {
            out << "[PASSED] NeighbourGrid cell size correct and optimal" << std::endl;
        }
    }

    bool Result::success() const {
        return voxelConflicts == 0 && biggerVoxelConflicts != 0 &&
               neighbourListConflicts == 0 && smallerNeighbourListConflicts != 0;
    }

    /* Generates random pair distant on x coordinate by distance and checks overlap */
    ShapePair random_pair(RND *rnd, double distance) {
        Vector<RSA_SPATIAL_DIMENSION> pos;
        pos[0] = distance / 2;

        auto shape1 = generate_randomly_oriented_shape(-pos, rnd);
        auto shape2 = generate_randomly_oriented_shape(pos, rnd);
        return ShapePair(std::move(shape1), std::move(shape2));
    }

    Result perform_test(unsigned long max_tries) {
        std::cout << "[INFO] Performing test for " << max_tries << " pairs" << std::endl;

        Result result;
        RND rnd;
        InfoLooper looper(max_tries, 10000, "pairs tested...");
        while (looper.step()) {
            const double SD_SQRT = std::sqrt(RSA_SPATIAL_DIMENSION);

            result.evaluateForVoxel(random_pair(&rnd, 0.9999*SD_SQRT * RSAShape::getVoxelSpatialSize()));
            result.evaluateForBiggerVoxel(random_pair(&rnd, 1.05*SD_SQRT * RSAShape::getVoxelSpatialSize()));
            result.evaluateForNeighbourList(random_pair(&rnd, 1.0001 * RSAShape::getNeighbourListCellSize()));
            result.evaluateForSmallerNeighbourList(random_pair(&rnd, 0.95 * RSAShape::getNeighbourListCellSize()));
        }
        return result;
    }
}

namespace shape_sizetest
{
    int main(int argc, char **argv) {
        if (argc < 5) {
            std::cerr << "Usage: ./rsa_test shape_sizetest [particle] [attibutes] [max_tries]" << std::endl;
            return TEST_ERROR;
        }

        ShapeFactory::initShapeClass(argv[2], argv[3]);
        unsigned long max_tries = std::stoul(argv[4]);
        if (max_tries <= 0) {
            std::cerr << "[ERROR] max_tries <= 0" << std::endl;
            return TEST_ERROR;
        }

        Result result = perform_test(max_tries);
        result.print(std::cout);
        return result.success() ? TEST_SUCCESS : TEST_FAILURE;
    }
}