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
#include "../rsa3d/ShapeFactory.h"
#include "../rsa3d/Utils.h"
#include "utility/MockBC.h"

namespace
{
    RSAShape *generateShape(RND *rnd, const Vector<RSA_SPATIAL_DIMENSION> &position) {
        auto shape = ShapeFactory::createShape(rnd);

        double arrayPos[RSA_SPATIAL_DIMENSION];
        position.copyToArray(arrayPos);
        shape->translate(arrayPos);

        double angularSize = RSAShape::getVoxelAngularSize();
        std::array<double, RSA_ANGULAR_DIMENSION> orientation{};
        std::for_each(orientation.begin(), orientation.end(), [rnd, angularSize](double &elem) {
            elem = rnd->nextValue() * angularSize;
        });
        shape->rotate(orientation);

        return shape;
    }

    ShapePairFactory::ShapePair generateShapePair(RND *rnd, double distance) {
        Vector<RSA_SPATIAL_DIMENSION> pos;
        pos[0] = distance / 2;
        return ShapePairFactory::ShapePair(generateShape(rnd, -pos), generateShape(rnd, pos));
    }
}

namespace shape_sizetest
{
    int main(int argc, char **argv) {
        if (argc < 5)   die("Usage: ./rsa_test shape_sizetest [particle] [attibutes] [max_tries]");

        ShapeFactory::initShapeClass(argv[2], argv[3]);
        unsigned long max_tries = std::stoul(argv[4]);
        if (max_tries <= 0)     die("Wrong input. Aborting.");

        std::size_t voxelConflicts{}, biggerVoxelConflicts{}, neighbourListConflicts{}, smallerNeighbourListConflicts{};
        RND rnd;
        for (std::size_t i = 0; i < max_tries; i++) {
            MockBC bc;
            auto voxelPair = generateShapePair(&rnd,
                    RSAShape::getVoxelSpatialSize() * 1.9999 * std::sqrt(RSA_SPATIAL_DIMENSION));
            auto neighbourPair = generateShapePair(&rnd, RSAShape::getNeighbourListCellSize() * 1.0001);
            auto biggerVoxelPair = generateShapePair(&rnd,
                    RSAShape::getVoxelSpatialSize() * 2.1 * std::sqrt(RSA_SPATIAL_DIMENSION));
            auto smallerNeighbourPair = generateShapePair(&rnd, RSAShape::getNeighbourListCellSize() * 0.95);

            if (!voxelPair.first()->overlap(&bc, voxelPair.second()))                   voxelConflicts++;
            if (neighbourPair.first()->overlap(&bc, neighbourPair.second()))            neighbourListConflicts++;
            if (!biggerVoxelPair.first()->overlap(&bc, voxelPair.second()))             biggerVoxelConflicts++;
            if (smallerNeighbourPair.first()->overlap(&bc, neighbourPair.second()))     smallerNeighbourListConflicts++;

            if ((i % 10000) == 9999)
                std::cout << (i + 1) << " pairs tested..." << std::endl;
        }

        std::cout << std::endl;
        std::cout << ">> Voxel overlap conflicts: " << voxelConflicts << std::endl;
        if (biggerVoxelConflicts == 0)
            std::cout << ">> Voxel size too small, no conflicts found in slightly bigger" << std::endl;
        else
            std::cout << ">> Voxel size optimal" << std::endl;

        std::cout << ">> NeighbourGrid cell overlap conflicts: " << neighbourListConflicts << std::endl;
        if (smallerNeighbourListConflicts == 0)
            std::cout << ">> NeighbourGrid cell too big, no conflicts found in slightly smaller" << std::endl;
        else
            std::cout << ">> NeighbourGrid cell size optimal" << std::endl;

        return EXIT_SUCCESS;
    }
}