//
// Created by PKua on 17.06.18.
//

#include "PackingOverlapsTest.h"
#include "utility/MockBC.h"
#include "../rsa3d/ShapeFactory.h"

bool pack_ovtest::test_packing_overlaps(const Packing *packing) {
    std::size_t size = packing->size();
    std::size_t total = size * (size - 1) / 2;
    std::cout << "[INFO] Checking " << size << " shapes, " << total << " checks required..." << std::endl;

    std::size_t counter = 0;
    std::size_t overlaps = 0;
    MockBC bc;

    for (std::size_t i = 0; i < size; i++) {
        for (std::size_t j = i + 1; j < size; j++) {
            counter++;
            auto shape1 = (*packing)[i];
            auto shape2 = (*packing)[j];
            if (shape1->overlap(&bc, shape2) != 0) {
                overlaps++;
            }

            std::size_t localCounter = counter;
            if (localCounter % 1000000 == 0) {
                std::cout << localCounter << " pairs tested..." << std::endl;
            }
        }
    }

    if (overlaps == 0)  std::cout << "[SUCCESS] ";
    else                std::cout << "[FAILURE] ";
    std::cout << "Finished. " << overlaps << " overlaps found" << std::endl;
    return overlaps == 0;
}

int pack_ovtest::main(int argc, char **argv) {
    if (argc < 4)   die("Usage: ./rsa_test pack_ovtest <config> <packing file>");

    Parameters params(argv[2]);
    ShapeFactory::initShapeClass(params.particleType, params.particleAttributes);
    std::string fileIn(argv[3]);

    auto packing = PackingGenerator::fromFile(fileIn);
    bool result = test_packing_overlaps(packing);
    delete packing;
    return result ? EXIT_SUCCESS : EXIT_FAILURE;
}
