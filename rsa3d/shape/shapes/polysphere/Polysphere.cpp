//
// Created by ciesla on 11/13/25.
//

#include "Polysphere.h"

std::vector<std::array<double, 3>> Polysphere::sphereCentre;
std::vector<double> Polysphere::sphereR;

/* If one wants to center automatically on the largest disk:
 * number_of_disks [xy|rt] c01 c02 r0 c11 c12 r1 c21 c22 r2 ... area
 *
 * If one wants to retain original centering:
 * number_of_disks [xy|rt] c01 c02 r0 c11 c12 r1 c21 c22 r2 ... area 'dontcenter' insphere_radius */
void Polysphere::initClass(const std::string &args){
    Polysphere::sphereCentre.clear();
    Polysphere::sphereR.clear();

    std::istringstream in(args);

    size_t n;
    in >> n;
    for (size_t i=0; i<n; i++) {
        double c1, c2, c3, sphereR_;
        in >> c1;
        in >> c2;
        in >> c3;
        in >> sphereR_;
        if (sphereR_ <= 0.0)
            throw std::runtime_error("sphereR <= 0 for " + std::to_string(i) + " coord");

        Polysphere::sphereCentre.push_back({c1, c2, c3});
        Polysphere::sphereR.push_back(sphereR_);
    }

    double insphereRadius;
    in >> insphereRadius;

    double initialArea;
    in >> initialArea;

    if (initialArea==0.0){
        initialArea = Polysphere::mcArea(100000000);
    }
    Polysphere::normalizeArea(initialArea);
    insphereRadius /= std::pow(initialArea,1.0/3.0);

    double circumsphereRadius = 0;
    for (size_t i=0; i<Polysphere::sphereCentre.size(); i++) {
        std::array<double, 3> coordinates = Polysphere::sphereCentre[i];
        double r = std::sqrt(
                coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] + coordinates[2] * coordinates[2]) +
                   Polysphere::sphereR[i];
        if (circumsphereRadius < r)
            circumsphereRadius = r;
    }
    Validate(insphereRadius <= circumsphereRadius);

    ShapeStaticInfo<3, 3> shapeInfo;

    shapeInfo.setCircumsphereRadius(circumsphereRadius);
    shapeInfo.setInsphereRadius(insphereRadius);
    shapeInfo.setAngularVoxelSize(2*M_PI);
    shapeInfo.setSupportsSaturation(true);
    shapeInfo.setDefaultCreateShapeImpl <Polysphere> ();

    Shape::setShapeStaticInfo(shapeInfo);
}

double Polysphere::mcArea(size_t mcTrials){
    double maxSize = 1.0;
    for (size_t i = 0; i < Polysphere::sphereCentre.size(); i++){
        std::array<double, 3> coordinates = Polysphere::sphereCentre[i];
        double r = std::sqrt(
                coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] + coordinates[2] * coordinates[2]) +
                   Polysphere::sphereR[i];
        if (maxSize < r)
            maxSize = r;
    }
    RND rnd(12345);
    double x, y, z, dx, dy, dz;
    size_t inside = 0;

    for (size_t i=0; i<mcTrials; i++){
        x = (2.0*rnd.nextValue()-1.0)*maxSize;
        y = (2.0*rnd.nextValue()-1.0)*maxSize;
        z = (2.0*rnd.nextValue()-1.0)*maxSize;
        for(size_t index = 0; index < Polysphere::sphereR.size(); index++){
            std::array<double, 3> coordinates = Polysphere::sphereCentre[i];
            dx = coordinates[0] - x;
            dy = coordinates[1] - y;
            dz = coordinates[2] - z;
            if (dx*dx + dy*dy +dz*dz < Polysphere::sphereR[index]*Polysphere::sphereR[index]){
                inside++;
                break;
            }
        }
    }
    return ((double)inside/(double)mcTrials)*8.0*maxSize*maxSize*maxSize;
}

void Polysphere::normalizeArea(double area){
    double denominator = std::pow(area,1.0/3.0);

    for(size_t index = 0; index < Polysphere::sphereR.size(); index++){
        for(int j=0; j<3; j++) {
            Polysphere::sphereCentre[index][j] /= denominator;
        }
        Polysphere::sphereR[index] /= denominator;
    }
}

std::array<double, 3> Polysphere::getStaticSpherePosition(std::size_t index, const Vector<3> position, const Orientation<3> orientation){
    Expects(index < Polysphere::sphereCentre.size());
    const double s_alpha = sin(orientation[0]);
    const double s_beta = sin(orientation[1]);
    const double s_gamma = sin(orientation[2]);
    const double c_alpha = cos(orientation[0]);
    const double c_beta = cos(orientation[1]);
    const double c_gamma = cos(orientation[2]);
    double x, y, z;
    x = position[0] + (
            Polysphere::sphereCentre[index][0] * (c_alpha * c_gamma - s_alpha * c_beta * s_gamma) +
            Polysphere::sphereCentre[index][1] * (-c_alpha * s_gamma - s_alpha * c_beta * c_gamma) +
            Polysphere::sphereCentre[index][2] * (s_alpha * s_beta)
    );
    y = position[1] + (
            Polysphere::sphereCentre[index][0] * (s_alpha * c_gamma + c_alpha * c_beta * s_gamma) +
            Polysphere::sphereCentre[index][1] * (-s_alpha * s_gamma + c_alpha * c_beta * c_gamma) +
            Polysphere::sphereCentre[index][2] * (-c_alpha * s_beta)
    );
    z = position[2] + (
            Polysphere::sphereCentre[index][0] * (s_beta * s_gamma) +
            Polysphere::sphereCentre[index][1] * (s_beta * c_gamma) +
            Polysphere::sphereCentre[index][2] * (c_beta)
    );
    return {x, y, z};
}

std::array<double, 3> Polysphere::getSpherePosition(size_t index) const {
    return Polysphere::getStaticSpherePosition(index, this->getPosition(), this->getOrientation());
}


bool Polysphere::sphereSphereIntersect(size_t sphere0, const Vector<3> shape0Position, const Orientation<3> shape0Orientation,
                                 size_t sphere1, const Vector<3> shape1Position, const Orientation<3> shape1Orientation){
    std::array<double, 3> coords0 = Polysphere::getStaticSpherePosition(sphere0, shape0Position, shape0Orientation);
    std::array<double, 3> coords1 = Polysphere::getStaticSpherePosition(sphere1, shape1Position, shape1Orientation);

    double distance2 =
            (coords0[0] - coords1[0])*(coords0[0] - coords1[0]) +
            (coords0[1] - coords1[1])*(coords0[1] - coords1[1]) +
            (coords0[2] - coords1[2])*(coords0[2] - coords1[2]);
    distance2 -= (Polysphere::sphereR[sphere0] + Polysphere::sphereR[sphere1])*(Polysphere::sphereR[sphere0] + Polysphere::sphereR[sphere1]);
    return distance2 < 0.0;
}

std::array<double, 2> Polysphere::minmaxCos(double theta, double dt){

    while (theta>2*M_PI)
        theta -= 2*M_PI;
    while (theta<0.0)
        theta += 2*M_PI;

    double v1 = std::cos(theta);
    double v2 = std::cos(theta+dt);
    std::array<double, 2> result{
                             std::min(v1, v2),
                             std::max(v1,v2)
                     };

    if ( (theta < 0 && theta+dt > 0 ) || (theta < 2.0*M_PI && theta+dt > 2.0*M_PI) ){
        result[1] = 1.0;
    }else if (theta < M_PI && theta+dt > M_PI){
        result[0] = -1.0;
    }
    return result;
}

std::array<double, 2> Polysphere::minmaxSin(double theta, double dt){

    while (theta>2*M_PI)
        theta -= 2*M_PI;
    while (theta<0.0)
        theta += 2*M_PI;

    double v1 = std::sin(theta);
    double v2 = std::sin(theta+dt);

    std::array<double, 2> result{
                             std::min(v1, v2),
                             std::max(v1,v2)
                     };
    if (theta < 0.5*M_PI && theta+dt > 0.5*M_PI ){
        result[1] = 1.0;
    }else if (theta < 1.5*M_PI && theta+dt > 1.5*M_PI){
        result[0] = -1.0;
    }
    return result;
}
// sphere0 is a part of a particle in a packing and sphere1 is a part of virtual particle placed anywhere in the voxel
// the method checks if a the sphere1 will intersect with sphere0 nevertheless the virtual particle is placed inside the voxel
bool Polysphere::sphereVoxelIntersect(size_t sphere0, const Vector<3> shape0Position, const Orientation<3> shape0Orientation,
                                  size_t sphere1, const Vector<3> voxelPosition, const Orientation<3> voxelOrientation,
                                  double spatialSize, double angularSize, std::array<std::array<double,2>, 6> minmaxTrigonometricArray){


    std::array<double, 2> sinAlpha = minmaxTrigonometricArray[0];
    std::array<double, 2> cosAlpha = minmaxTrigonometricArray[1];
    std::array<double, 2> sinBeta = minmaxTrigonometricArray[2];
    std::array<double, 2> cosBeta = minmaxTrigonometricArray[3];
    std::array<double, 2> sinGamma = minmaxTrigonometricArray[4];
    std::array<double, 2> cosGamma = minmaxTrigonometricArray[5];

    double xmin = voxelPosition[0] +
            (
            Polysphere::sphereCentre[sphere1][0] * (cosAlpha[0] * cosGamma[0] - sinAlpha[1] * cosBeta[1] * sinGamma[1]) +
            Polysphere::sphereCentre[sphere1][1] * (-cosAlpha[1] * sinGamma[1] - sinAlpha[1] * cosBeta[1] * cosGamma[1]) +
            Polysphere::sphereCentre[sphere1][2] * (sinAlpha[0] * sinBeta[0])
            );
    double xmax = voxelPosition[0] + spatialSize +
            (
            Polysphere::sphereCentre[sphere1][0] * (cosAlpha[1] * cosGamma[1] - sinAlpha[0] * cosBeta[0] * sinGamma[0]) +
            Polysphere::sphereCentre[sphere1][1] * (-cosAlpha[0] * sinGamma[0] - sinAlpha[0] * cosBeta[0] * cosGamma[0]) +
            Polysphere::sphereCentre[sphere1][2] * (sinAlpha[1] * sinBeta[1])
            );

    double ymin = voxelPosition[1] +
            (
            Polysphere::sphereCentre[sphere1][0] * (sinAlpha[0] * cosGamma[0] + cosAlpha[0] * cosBeta[0] * sinGamma[0]) +
            Polysphere::sphereCentre[sphere1][1] * (- sinAlpha[1] * sinGamma[1] + cosAlpha[0] * cosBeta[0] * cosGamma[0]) +
            Polysphere::sphereCentre[sphere1][2] * (- cosAlpha[1] * sinBeta[1])
            );
    double ymax = voxelPosition[1] + spatialSize +
            (
            Polysphere::sphereCentre[sphere1][0] * (sinAlpha[1] * cosGamma[1] + cosAlpha[1] * cosBeta[1] * sinGamma[1]) +
            Polysphere::sphereCentre[sphere1][1] * (- sinAlpha[0] * sinGamma[0] + cosAlpha[1] * cosBeta[1] * cosGamma[1]) +
            Polysphere::sphereCentre[sphere1][2] * (- cosAlpha[0] * sinBeta[0])
            );

    double zmin = voxelPosition[2] +
            (
            Polysphere::sphereCentre[sphere1][0] * (sinBeta[0] * sinGamma[0]) +
            Polysphere::sphereCentre[sphere1][1] * (sinBeta[0] * cosGamma[0]) +
            Polysphere::sphereCentre[sphere1][2] * (cosBeta[0])
            );
    double zmax = voxelPosition[2] + spatialSize +
            (
            Polysphere::sphereCentre[sphere1][0] * (sinBeta[1] * sinGamma[1]) +
            Polysphere::sphereCentre[sphere1][1] * (sinBeta[1] * cosGamma[1]) +
            Polysphere::sphereCentre[sphere1][2] * (cosBeta[1])
            );

    std::array<double, 3> coordinates = Polysphere::getStaticSpherePosition(sphere0, shape0Position, shape0Orientation);

    double xmax2 = std::max( (xmin - coordinates[0])*(xmin - coordinates[0]) ,(xmax - coordinates[0])*(xmax - coordinates[0]));
    double ymax2 = std::max( (ymin - coordinates[1])*(ymin - coordinates[1]) ,(ymax - coordinates[1])*(ymax - coordinates[1]));
    double zmax2 = std::max( (zmin - coordinates[2])*(zmin - coordinates[2]) ,(zmax - coordinates[2])*(zmax - coordinates[2]));
    double r2 = (Polysphere::sphereR[sphere0] + Polysphere::sphereR[sphere1]) * (Polysphere::sphereR[sphere0] + Polysphere::sphereR[sphere1]);
    if (r2 > xmax2 + ymax2 + zmax2)
        return true;
    else
        return false;
}


double Polysphere::getVolume(unsigned short dim) const {
    if (dim != 3)
        throw std::runtime_error ("Polysphere supports only 3D packings");

    return 1.0;
}

bool Polysphere::overlap(BoundaryConditions<3> *bc, const Shape<3, 3> *s) const{
    switch (this->overlapEarlyRejection(bc, s)) {
        case TRUE:      return true;
        case FALSE:     return false;
        case UNKNOWN:   break;
    }

    Polysphere pol = dynamic_cast<const Polysphere&>(*s);
    this->applyBC(bc, &pol);

    Vector<3> polPosition = pol.getPosition();
    Orientation<3> polOrientation = pol.getOrientation();

    Vector<3> thisPosition = this->getPosition();
    Orientation<3> thisOrientation = this->getOrientation();

    for (size_t i = 0; i < Polysphere::sphereCentre.size(); i++){
        for (size_t j = 0; j < Polysphere::sphereCentre.size(); j++){
            if ( Polysphere::sphereSphereIntersect(i, polPosition, polOrientation, j, thisPosition, thisOrientation) )
                return true;
        }
    }
    return false;
}

bool Polysphere::fullAngleVoxelInside(BoundaryConditions<3> *bc, const Vector<3> &voxelPosition,
                                    double spatialSize) const
{
    double voxelCircumsphereRadius = spatialSize / M_SQRT2;
    double insphereRadius = Polysphere::getInsphereRadius();

    for (std::size_t i = 0; i < Polysphere::sphereR.size(); i++) {
        double voxelDistance2 = bc->distance2(this->getSpherePosition(i), voxelPosition);
        if (voxelDistance2 <= std::pow(Polysphere::sphereR[i] + insphereRadius - voxelCircumsphereRadius, 2))
            return true;
    }
    return false;
}

bool Polysphere::voxelInside(BoundaryConditions<3> *bc, const Vector<3> &voxelPosition, const Orientation<3> &voxelOrientation, double spatialSize, double angularSize) const{

    double angularVoxelSize = Shape<3, 3>::getAngularVoxelSize();
    if (
            voxelOrientation[0] > angularVoxelSize ||
            voxelOrientation[1] > angularVoxelSize ||
            voxelOrientation[2] > angularVoxelSize
    )
        return true;

    switch(this->voxelInsideEarlyRejection(bc, voxelPosition, voxelOrientation, spatialSize, angularSize)) {
        case TRUE:      return true;
        case FALSE:     return false;
        case UNKNOWN:   break;
    }

    // Version optimized for full angle - it was made mainly for OrientedFibrinogen, because it took ages to generate
    if (angularSize >= 2*M_PI)
        return this->fullAngleVoxelInside(bc, voxelPosition, spatialSize);

    Vector<3> translation = bc->getTranslation(voxelPosition, this->getPosition());
    Vector<3> thisPosition = this->getPosition() + translation;
    Orientation<3> thisOrientation = this->getOrientation();

    std::array<std::array<double, 2>, 6> minMaxTrigonometricArray = {
            Polysphere::minmaxSin(voxelOrientation[0], angularSize),
            Polysphere::minmaxCos(voxelOrientation[0], angularSize),
            Polysphere::minmaxSin(voxelOrientation[1], angularSize),
            Polysphere::minmaxCos(voxelOrientation[1], angularSize),
            Polysphere::minmaxSin(voxelOrientation[2], angularSize),
            Polysphere::minmaxCos(voxelOrientation[2], angularSize)
    };

    // loop over disks in this particle
    for (size_t i = 0; i < Polysphere::sphereCentre.size(); i++){
        // loop over disks in virtual particle inside voxel
        for (size_t j = 0; j < Polysphere::sphereCentre.size(); j++){
            // if a disk of virtual particle (placed anywhere in the voxel) intersects with disk of this particle, the voxel is fully inside the exclusion zone
            if (Polysphere::sphereVoxelIntersect(i, thisPosition, thisOrientation, j, voxelPosition, voxelOrientation, spatialSize, angularSize, minMaxTrigonometricArray) )
                return true;
        }
    }
    return false;
}

Shape<3, 3> *Polysphere::clone() const {
    return new Polysphere(*this);
}

std::string Polysphere::toPovray() const{
    std::stringstream out;
    out.precision(std::numeric_limits< double >::max_digits10);
    for (size_t i=0; i < Polysphere::sphereCentre.size(); i++) {
        std::array<double, 3> spherePosition = this->getSpherePosition(i);
        double r = Polysphere::sphereR[i];
        out << "  sphere { < " << spherePosition[0] << ", " << spherePosition[1] << ", " << spherePosition[2] << ">, ";
        out << r << std::endl;
        out << "    texture { finish { ambient 1 diffuse 0 } pigment { color Red} }" << std::endl << "}" << std::endl;
    }

    return out.str();
}

std::string Polysphere::toWolfram() const {
    std::stringstream out;
    out.precision(std::numeric_limits<double>::max_digits10);

    out << std::fixed << "{";
    for (std::size_t i = 0, max = Polysphere::sphereCentre.size(); i < max; i++) {
        std::array<double, 3> spherePosition = this->getSpherePosition(i);
        out << "Sphere[{" << spherePosition[0] << ", " << spherePosition[1] << ", " << spherePosition[2] << "}, " << Polysphere::sphereR[i] << "]";
        if (i < max - 1)
            out << ", ";
    }
    out << "}";

    return out.str();
}