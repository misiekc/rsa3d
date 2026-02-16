//
// Created by ciesla on 15.11.2025.
//

#include "Kmer3D.h"

void Kmer3D::initClass(const std::string &attr) {
    std::istringstream attrStream(attr);
    int numberOfSpheres;
    double distance;
    attrStream >> numberOfSpheres >> distance;

    ValidateMsg(attrStream, "arguments format: [number of disks] [distance between neighbouring disks]");
    Validate(numberOfSpheres >= 2);
    Validate(distance >= 0);
    Validate(distance < 2*std::sqrt(2));
    Validate(distance <= 2*M_SQRT2);

    Polysphere::initClass(Kmer3D::preparePolysphereAttr(numberOfSpheres, distance));
    ShapeStaticInfo<3,3> shapeInfo = Shape::getShapeStaticInfo();
    shapeInfo.setDefaultCreateShapeImpl <Kmer3D> ();
//    shapeInfo.setAngularVoxelSize({2*M_PI, 1.0, 2*M_PI});
    Shape::setShapeStaticInfo(shapeInfo);
}

std::string Kmer3D::preparePolysphereAttr(int numberOfDisks, double distance) {
    double volume = Kmer3D::calculateVolume(numberOfDisks, distance);
    double inSphereRadius = ((numberOfDisks % 2)==0)? 0.5*std::sqrt(4-0.5*distance*distance) : 1.0;

    std::ostringstream polydiskAttrStream;
    polydiskAttrStream << numberOfDisks;
    double firstDiskOffset = distance * (numberOfDisks - 1) / 2;
    for (int i = 0; i < numberOfDisks; i++)
        polydiskAttrStream << " " << (distance * i - firstDiskOffset) << " 0 0 1 ";
    polydiskAttrStream << inSphereRadius << " " << volume;

    return polydiskAttrStream.str();
}

double Kmer3D::calculateVolume(int numberOfSpheres, double distance) {
    // length of particle
    double length = (numberOfSpheres-1)*distance + 2.0;
    // volume of a dimer containing two spheres at a distance x (between their centers)
    double vdim = (length < 4.0) ? M_PI * (12.0 * distance - distance * distance * distance + 16.0) / 12.0 : (8.0 / 3.0) * M_PI;
    // volume between spherocylinder and two neighbouring spheres
    double vrest = ((4.0 / 3.0) + (length-2.0)) * M_PI - vdim;

    if (distance < 2.0)
        // volume of spherocylinder - (numberOfSpheres-1)*vrest
        return M_PI * (length-2.0 + (4.0 / 3.0)) - (numberOfSpheres - 1) * vrest;
    else
        return (4.0 / 3.0) * M_PI * numberOfSpheres;
}

std::array<Matrix<3,3>, 2> Kmer3D::getMinMaxMatrices(const Orientation<3> &voxelOrientation, const Orientation<3> &angularSize) const{
    Orientation<3> translatedVoxelOrientation = Polysphere::fromVoxel(voxelOrientation);
    double angularSize1 = std::asin(voxelOrientation[1]+angularSize[1]-1.0) - std::asin(voxelOrientation[1]-1.0);

    Orientation<3> sinTheta{ std::sin(translatedVoxelOrientation[0]), std::sin(translatedVoxelOrientation[1]), std::sin(translatedVoxelOrientation[2]) };
    Orientation<3> cosTheta{ std::cos(translatedVoxelOrientation[0]), std::cos(translatedVoxelOrientation[1]), std::cos(translatedVoxelOrientation[2]) };
    Orientation<3> sinDt{ std::sin(angularSize[0]), std::sin(angularSize1), std::sin(angularSize[2])};
    Orientation<3> cosDt{ std::cos(angularSize[0]), std::cos(angularSize1), std::cos(angularSize[2])};


    std::array<std::array<double, 2>, 2> minmaxSinCos = Polysphere::minmaxSinCos(translatedVoxelOrientation[1], angularSize1, sinTheta[1], cosTheta[1], sinDt[1], cosDt[1]);
    std::array<double, 2> sin_ay = minmaxSinCos[0];
    std::array<double, 2> cos_ay = minmaxSinCos[1];
    minmaxSinCos = Polysphere::minmaxSinCos(translatedVoxelOrientation[2], angularSize[2], sinTheta[2], cosTheta[2], sinDt[2], cosDt[2]);
    std::array<double, 2> sin_az = minmaxSinCos[0];
    std::array<double, 2> cos_az = minmaxSinCos[1];

    Matrix<3, 3, double> minMatrix({
        cos_ay[0] * cos_az[0],
        -sin_az[1],
        sin_ay[0] * cos_az[0],

        cos_ay[0] * sin_az[0],
        cos_az[0],
        sin_ay[0] * sin_az[0],

        -sin_ay[1],
        0.0,
        cos_ay[0]
    });
    Matrix<3, 3, double> maxMatrix({
        cos_ay[1] * cos_az[1],
        -sin_az[0],
        sin_ay[1] * cos_az[1],

        cos_ay[1] * sin_az[1],
        cos_az[1],
        sin_ay[1] * sin_az[1],

        -sin_ay[0],
        0.0,
        cos_ay[1]
    });
    return{minMatrix, maxMatrix};
}

bool Kmer3D::voxelInside(BoundaryConditions<3> *bc, const Vector<3> &voxelPosition, const Orientation<3> &voxelOrientation, double spatialSize, const Orientation<3> &angularSize) const{
//    if (voxelOrientation[0] > 0.0)
    if (voxelOrientation[0] > 0.0 || voxelOrientation[1]>1.0)
        return true;
    else
        return Polysphere::voxelInside(bc, voxelPosition, voxelOrientation, spatialSize, angularSize);
}
