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
#include "../rsa3d/Utils.h"
#include "utility/MockBC.h"
#include "utility/InfoLooper.h"
#include "utility/ShapeGenerators.h"

namespace
{
    struct Result {
        std::size_t voxelConflicts{};
        std::string voxelConflictExample = "";

        std::size_t biggerVoxelConflicts{};

        std::size_t neighbourListConflicts{};
        std::string neighbourListConflictExample = "";

        std::size_t smallerNeighbourListConflicts{};

        
        void print(std::ostream &out) const;
        bool success() const;
    };

    /* Presents result on the out ostream */
    void Result::print(std::ostream &out) const {
        out << std::endl;
        if (voxelConflicts != 0){
            out << "[FAILED] Voxel size too big. Empty spaces will be produced" << std::endl;
            out << "Example configuration of two non-overlapping (according to the overlap() method) shapes:" << std::endl << voxelConflictExample << std::endl << std::endl;
        }
        else if (biggerVoxelConflicts == 0){
            out << "[WARNING] Voxel size too small, no conflicts found in a slightly bigger. Check more pairs before enlarging" << std::endl;
        }
        else
            out << "[PASSED] Voxel size correct and optimal" << std::endl;

        if (neighbourListConflicts != 0){
            out << "[FAILED] NeighbourGrid cell to small. Overlapping shapes will occur" << std::endl;
            out << "Example configuration of two overlapping (according to the overlap() method) shapes:" << std::endl << neighbourListConflictExample << std::endl << std::endl;
        }
        else if (smallerNeighbourListConflicts == 0){
            out << "[WARNING] NeighbourGrid cell too big, no conflicts found in a slightly smaller. Check more pairs before reducing" << std::endl;
        }
        else
            out << "[PASSED] NeighbourGrid cell size correct and optimal" << std::endl;
    }

    bool Result::success() const {
        return neighbourListConflicts == 0;
    }

    std::unique_ptr<RSAShape> shape1, shape2;

    /* Generates random pair distant on x coordinate by distance and checks overlap */
    bool random_pair_overlap(RND *rnd, double distance) {
        Vector<RSA_SPATIAL_DIMENSION> pos;
        pos[0] = distance / 2;
        shape1 = generate_randomly_oriented_shape(-pos, rnd);
        shape2 = generate_randomly_oriented_shape(pos, rnd);

        RSAMockBC bc;
        return shape1->overlap(&bc, shape2.get());
    }

    Result perform_test(unsigned long max_tries) {
        std::cout << "[INFO] Performing test for " << max_tries << " pairs" << std::endl;

        Result result;
        RND rnd;
        InfoLooper looper(max_tries, 10000, "pairs tested...");
        while (looper.step()) {
            if (!random_pair_overlap(&rnd, RSAShape::getVoxelSpatialSize() * 0.9999 * std::sqrt(RSA_SPATIAL_DIMENSION))){
                result.voxelConflicts++;
                if (result.voxelConflictExample.empty()){
                	result.voxelConflictExample = "Graphics3D[{ \r\n" +
                								   shape1->toWolfram() + ", \r\n" +
												   shape2->toWolfram() +
												   "\r\n }]" ;
                }
            }
            if (!random_pair_overlap(&rnd, RSAShape::getVoxelSpatialSize() * 1.05 * std::sqrt(RSA_SPATIAL_DIMENSION))){
                result.biggerVoxelConflicts++;
            }
            if (random_pair_overlap(&rnd, RSAShape::getNeighbourListCellSize() * 1.0001)){
                result.neighbourListConflicts++;
                if (result.neighbourListConflictExample.empty()){
                	result.neighbourListConflictExample = "Graphics3D[{ \r\n" +
                								   shape1->toWolfram() + ", \r\n" +
												   shape2->toWolfram() +
												   "\r\n }]" ;
                }
            }
            if (random_pair_overlap(&rnd, RSAShape::getNeighbourListCellSize() * 0.95)){
                result.smallerNeighbourListConflicts++;
            }
        }
        return result;
    }
}

namespace shape_sizetest
{
    int main(int argc, char **argv) {
        if (argc < 5)   die("Usage: ./rsa_test shape_sizetest [particle] [attibutes] [max_tries]");

        ShapeFactory::initShapeClass(argv[2], argv[3]);
        unsigned long max_tries = std::stoul(argv[4]);
        if (max_tries <= 0)     die("Wrong input. Aborting.");

        Result result = perform_test(max_tries);
        result.print(std::cout);
        return result.success() ? EXIT_SUCCESS : EXIT_FAILURE;
    }
}
