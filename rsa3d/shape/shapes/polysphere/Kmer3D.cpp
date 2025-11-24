//
// Created by pkua on 15.11.2019.
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
    Validate(distance <= 2*M_SQRT2);

    Polysphere::initClass(Kmer3D::preparePolysphereAttr(numberOfSpheres, distance));
    ShapeStaticInfo<3,3> shapeInfo = Shape::getShapeStaticInfo();
    shapeInfo.setDefaultCreateShapeImpl <Kmer3D> ();
    Shape::setShapeStaticInfo(shapeInfo);
}

std::string Kmer3D::preparePolysphereAttr(int numberOfDisks, double length) {
    double volume = Kmer3D::calculateVolume(numberOfDisks, length);

    std::ostringstream polydiskAttrStream;
    polydiskAttrStream << numberOfDisks;
    double firstDiskOffset = length * (numberOfDisks - 1) / 2;
    for (int i = 0; i < numberOfDisks; i++)
        polydiskAttrStream << " 0 0 " << (length * i - firstDiskOffset) << " 1";
    polydiskAttrStream << " 1 " << volume;

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

std::array<std::array<double, 2>, 3> Kmer3D::getMinMaxVoxelCoordinates(size_t sphereIndex, const Vector<3> &position, const Orientation<3> &orientation, double spatialSize, double angularSize) const{
    std::array<double, 2>   sinAlpha = Polysphere::minmaxSin(orientation[0], angularSize),
                            cosAlpha = Polysphere::minmaxCos(orientation[0], angularSize),
                            sinBeta  = Polysphere::minmaxSin(orientation[1], angularSize),
                            cosBeta  = Polysphere::minmaxCos(orientation[1], angularSize);

    double xmin = position[0] +
        (
            Polysphere::sphereCentre[sphereIndex][0] * (cosAlpha[0]) +
            Polysphere::sphereCentre[sphereIndex][1] * (-sinAlpha[1] * cosBeta[1]) +
            Polysphere::sphereCentre[sphereIndex][2] * (sinAlpha[0] * sinBeta[0])
        );
    double xmax = position[0] + spatialSize +
        (
            Polysphere::sphereCentre[sphereIndex][0] * (cosAlpha[1]) +
            Polysphere::sphereCentre[sphereIndex][1] * (-sinAlpha[0] * cosBeta[0]) +
            Polysphere::sphereCentre[sphereIndex][2] * (sinAlpha[1] * sinBeta[1])
        );

    double ymin = position[1] +
        (
            Polysphere::sphereCentre[sphereIndex][0] * (sinAlpha[0]) +
            Polysphere::sphereCentre[sphereIndex][1] * (cosAlpha[0] * cosBeta[0]) +
            Polysphere::sphereCentre[sphereIndex][2] * (- cosAlpha[1] * sinBeta[1])
        );

    double ymax = position[1] + spatialSize +
        (
            Polysphere::sphereCentre[sphereIndex][0] * (sinAlpha[1]) +
            Polysphere::sphereCentre[sphereIndex][1] * (cosAlpha[1] * cosBeta[1]) +
            Polysphere::sphereCentre[sphereIndex][2] * (- cosAlpha[0] * sinBeta[0])
        );

    double zmin = position[2] +
        (
            Polysphere::sphereCentre[sphereIndex][1] * (sinBeta[0]) +
            Polysphere::sphereCentre[sphereIndex][2] * (cosBeta[0])
        );

    double zmax = position[2] + spatialSize +
        (
            Polysphere::sphereCentre[sphereIndex][1] * (sinBeta[1]) +
            Polysphere::sphereCentre[sphereIndex][2] * (cosBeta[1])
        );

    return std::array<std::array<double, 2>, 3>{std::array<double, 2>{xmin, xmax}, std::array<double, 2>{ymin, ymax}, std::array<double, 2>{zmin, zmax}};;
}

bool Kmer3D::voxelInside(BoundaryConditions<3> *bc, const Vector<3> &voxelPosition, const Orientation<3> &voxelOrientation, double spatialSize, double angularSize) const{
    if (voxelOrientation[2] > 0.0)
        return true;
    else
        return Polysphere::voxelInside(bc, voxelPosition, voxelOrientation, spatialSize, angularSize);
}
