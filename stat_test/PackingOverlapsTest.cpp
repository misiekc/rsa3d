//
// Created by PKua on 17.06.18.
//

#include "PackingOverlapsTest.h"
#include "../rsa3d/FreeBC.h"
#include "../rsa3d/shape/ShapeFactory.h"
#include "../rsa3d/Parameters.h"
#include "../rsa3d/utils/OMPMacros.h"
#include "utils/ParalellInfoLooper.h"


bool pack_ovtest::test_packing_overlaps(const Packing &packing) {
    std::size_t size = packing.size();
    std::size_t total = size * (size - 1) / 2;
    std::cout << "[INFO] Checking " << size << " shapes, " << total << " checks required..." << std::endl;

    std::size_t overlaps = 0;
    RSAFreeBC bc;
    ParallelInfoLooper infoLooper(10000000);

    _OMP_PARALLEL_FOR
    for (std::size_t i = 0; i < size; i++) {
        for (std::size_t j = i + 1; j < size; j++) {
            auto shape1 = packing[i];
            auto shape2 = packing[j];
            if (shape1->overlap(&bc, shape2) != 0) {
                _OMP_ATOMIC
                overlaps++;
            }
            infoLooper.step();
        }
    }

    if (overlaps == 0)  std::cout << "[SUCCESS] ";
    else                std::cout << "[FAILURE] ";
    std::cout << "Finished. " << overlaps << " overlaps found" << std::endl;
    return overlaps != 0;
}

int pack_ovtest::main(int argc, char **argv) {
    if (argc < 4)   die("Usage: ./rsa_test pack_ovtest [config file] [packing file]");

    Parameters params(argv[2]);
    ShapeFactory::initShapeClass(params.particleType, params.particleAttributes);
    std::string fileIn(argv[3]);

    Packing packing;
    packing.restore(fileIn);
    bool overlapsFound = test_packing_overlaps(packing);
    return overlapsFound ? EXIT_FAILURE : EXIT_SUCCESS;
}
