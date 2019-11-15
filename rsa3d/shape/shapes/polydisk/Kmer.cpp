//
// Created by pkua on 15.11.2019.
//

#include "Kmer.h"

void Kmer::initClass(const std::string &attr) {
    std::istringstream attrStream(attr);
    int numberOfDisks;
    double distance;
    attrStream >> numberOfDisks >> distance;

    ValidateMsg(attrStream, "arguments format: [number of disks] [distance between disks]");
    Validate(numberOfDisks >= 2);
    Validate(distance >= 0);
    Validate(distance <= 2*M_SQRT2);

    Polydisk::initClass(Kmer::preparePolydiskAttr(numberOfDisks, distance));
}

std::string Kmer::preparePolydiskAttr(int numberOfDisks, double distance) {
    double area = calculateArea(numberOfDisks, distance);

    std::ostringstream polydiskAttrStream;
    polydiskAttrStream << numberOfDisks << " xy 0 0 1";
    for (int i = 1; i < numberOfDisks; i++)
        polydiskAttrStream << " " << (distance * i) << " 0 1";
    polydiskAttrStream << " " << area;

    return polydiskAttrStream.str();
}

double Kmer::calculateArea(int numberOfDisks, double distance) {
    Assert(distance >= 0);
    Assert(distance <= 2 * M_SQRT2);

    // Area based on Cieśla, M., Paja̧k, G., & Ziff, R. M., J. Chem. Phys., 2016, 145(4), 044708
    double S_pol{};
    double S_ad{};
    if (distance < 2) {
        S_pol = M_PI + (numberOfDisks - 1)*(distance/2*sqrt(4 - distance*distance) + 2*asin(distance/2));
        S_ad = distance/4*(std::sqrt(16 - distance*distance) - std::sqrt(4 - distance*distance)) - asin(distance/2);
    } else {
        S_pol = numberOfDisks*M_PI;
        S_ad = distance/4*(sqrt(16 - distance*distance)) - M_PI/2;
    }

    return S_pol + 2*(numberOfDisks - 1)*S_ad;
}
