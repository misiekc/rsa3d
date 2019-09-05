//
// Created by pkua on 04.09.2019.
//

#include <sstream>

#include "Fibrinogen.h"

/**
 * @brief Initializes Polydisk class so it represents linear fibrinogen particle, consisting in medium disk connected
 * with 2 large disks at the ends by chains of small disks.
 * @param attr parameters in format: (radius of large disks) (radius of a medium disk) (radius of small disks)
 * (number of small disks in one chain)
 */
void Fibrinogen::initClass(const std::string &attr) {
    std::istringstream attrStream;
    if (attr.empty())
        attrStream = std::istringstream("6.7 5.3 1.5 10");
    else
        attrStream = std::istringstream(attr);

    double largeDiskRadius{}, mediumDiskRadius{}, smallDiskRadius{};
    std::size_t numberOfSmallDisks{};
    attrStream >> largeDiskRadius >> mediumDiskRadius >> smallDiskRadius >> numberOfSmallDisks;
    ValidateMsg(attrStream, "Wrong format, expecting: (radius of large disks) (radius of a medium disk) "
                            "(radius of small disks) (radius of small disks in one chain)");
    Validate(largeDiskRadius > 0);
    Validate(mediumDiskRadius > 0);
    Validate(smallDiskRadius > 0);
    Validate(numberOfSmallDisks > 0);

    /* Polydisk attr format: number_of_disks [xy|rt] c01 c02 r0 c11 c12 r1 c21 c22 r2 ... area */
    std::ostringstream polydiskAttrStream;
    std::size_t numberOfDisks = 3 + 2*numberOfSmallDisks;
    polydiskAttrStream << numberOfDisks << " xy 0 0 " << mediumDiskRadius;

    double lastCircleDisplacement = mediumDiskRadius;
    for (std::size_t i = 0; i < numberOfSmallDisks; i++) {
        lastCircleDisplacement += smallDiskRadius;
        polydiskAttrStream << " " << lastCircleDisplacement << " 0 " << smallDiskRadius;
        polydiskAttrStream << " -" << lastCircleDisplacement << " 0 " << smallDiskRadius;
        lastCircleDisplacement += smallDiskRadius;
    }

    lastCircleDisplacement += largeDiskRadius;
    polydiskAttrStream << " " << lastCircleDisplacement << " 0 " << largeDiskRadius;
    polydiskAttrStream << " -" << lastCircleDisplacement << " 0 " << largeDiskRadius;

    double area = M_PI * (2*std::pow(largeDiskRadius, 2) + std::pow(mediumDiskRadius, 2)
                          + 2*numberOfSmallDisks * std::pow(smallDiskRadius, 2));
    polydiskAttrStream << " " << area;

    Polydisk::initClass(polydiskAttrStream.str());
}
