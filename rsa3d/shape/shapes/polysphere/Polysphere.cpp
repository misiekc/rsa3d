//
// Created by ciesla on 11/13/25.
//

#include "Polysphere.h"

std::vector<Vector<3>> Polysphere::sphereCentre;
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

        Polysphere::sphereCentre.push_back(Vector<3>({c1, c2, c3}));
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
        Vector<3> coordinates = Polysphere::sphereCentre[i];
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
    shapeInfo.setAngularVoxelSize({2*M_PI, 2, 2*M_PI});
    shapeInfo.setSupportsSaturation(true);
    shapeInfo.setDefaultCreateShapeImpl <Polysphere> ();

    Shape::setShapeStaticInfo(shapeInfo);
}

double Polysphere::mcArea(size_t mcTrials){
    double maxSize = 1.0;
    for (size_t i = 0; i < Polysphere::sphereCentre.size(); i++){
        Vector<3> coordinates = Polysphere::sphereCentre[i];
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
            Vector<3> coordinates = Polysphere::sphereCentre[i];
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

Vector<3> Polysphere::getStaticSpherePosition(std::size_t index, const Vector<3> &position, const Orientation<3> &orientation) {
    Matrix<3,3> rotation = Matrix<3,3>::rotation(orientation[0], orientation[1], orientation[2]);
    Vector<3> result = position + rotation * Polysphere::sphereCentre[index];
    return result;
}

/*
Vector<3> Polysphere::getStaticSpherePosition(std::size_t index, const Vector<3> &position, const Orientation<3> &orientation){
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
    return Vector<3>({x, y, z});
}
*/

Vector<3> Polysphere::getSpherePosition(size_t index) const {
    return Polysphere::getStaticSpherePosition(index, this->getPosition(), this->getOrientation());
}


bool Polysphere::sphereSphereIntersect(size_t sphere0, const Vector<3> &shape0Position, const Orientation<3> &shape0Orientation,
                                 size_t sphere1, const Vector<3> &shape1Position, const Orientation<3> &shape1Orientation){
    Vector<3> coords0 = Polysphere::getStaticSpherePosition(sphere0, shape0Position, shape0Orientation);
    Vector<3> coords1 = Polysphere::getStaticSpherePosition(sphere1, shape1Position, shape1Orientation);

    double distance2 = (coords0 - coords1).norm2();
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

Orientation<3> Polysphere::fromVoxel(const Orientation<3> &voxelOrientation) {
    Orientation<3> result;
    result[0] = voxelOrientation[0];
    result[1] = std::asin(voxelOrientation[1]-1.0);
    result[2] = voxelOrientation[2];
    return result;
}

double Polysphere::getVolume(unsigned short dim) const {
    if (dim != 3)
        throw std::runtime_error ("Polysphere supports only 3D packings");

    return 1.0;
}

void Polysphere::rotate(const Orientation<3> &voxelOrientation){
    Orientation<3> angles = this->getOrientation();
    Orientation<3> translatefOrientation = Polysphere::fromVoxel(voxelOrientation);
    angles[0] += translatefOrientation[0];
    angles[1] += translatefOrientation[1];
    angles[2] += translatefOrientation[2];
    this->setOrientation(angles);
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
std::array<Vector<3>, 2> Polysphere::getMinMaxVoxelCoordinates(size_t sphereIndex, const Vector<3> &voxelPosition, const Orientation<3> &voxelOrientation, double spatialSize, const Orientation<3> &angularSize) const {
    Orientation<3> orientation = Polysphere::fromVoxel(voxelOrientation);
    double angularSize1 = std::asin(voxelOrientation[1]+angularSize[1]-1.0) - std::asin(voxelOrientation[1]-1.0);

    std::array<double, 2> sin_ax = Polysphere::minmaxSin(orientation[0], angularSize[0]);
    std::array<double, 2> sin_ay = Polysphere::minmaxSin(orientation[1], angularSize1);
    std::array<double, 2> sin_az = Polysphere::minmaxSin(orientation[2], angularSize[2]);
    std::array<double, 2> cos_ax = Polysphere::minmaxCos(orientation[0], angularSize[0]);
    std::array<double, 2> cos_ay = Polysphere::minmaxCos(orientation[1], angularSize1);
    std::array<double, 2> cos_az = Polysphere::minmaxCos(orientation[2], angularSize[2]);

    Matrix<3, 3, double> minMatrix({
        cos_ay[0] * cos_az[0],
        sin_ax[0] * sin_ay[0] * cos_az[0] - cos_ax[1] * sin_az[1],
        cos_ax[0] * sin_ay[0] * cos_az[0] + sin_ax[0] * sin_az[0],

        cos_ay[0] * sin_az[0],
        sin_ax[0] * sin_ay[0] * sin_az[0] + cos_ax[0] * cos_az[0],
        cos_ax[0] * sin_ay[0] * sin_az[0] - sin_ax[1] * cos_az[1],

        -sin_ay[1],
        sin_ax[0] * cos_ay[0],
        cos_ax[0] * cos_ay[0]
    });
    Matrix<3, 3, double> maxMatrix({
        cos_ay[1] * cos_az[1],
        sin_ax[1] * sin_ay[1] * cos_az[1] - cos_ax[0] * sin_az[0],
        cos_ax[1] * sin_ay[1] * cos_az[1] + sin_ax[1] * sin_az[1],

        cos_ay[1] * sin_az[1],
        sin_ax[1] * sin_ay[1] * sin_az[1] + cos_ax[1] * cos_az[1],
        cos_ax[1] * sin_ay[1] * sin_az[1] - sin_ax[0] * cos_az[0],

        -sin_ay[0],
        sin_ax[1] * cos_ay[1],
        cos_ax[1] * cos_ay[1]
    });
    Vector<3> vMin = voxelPosition + minMatrix*Polysphere::sphereCentre[sphereIndex];
    Vector<3> vMax = voxelPosition + maxMatrix*Polysphere::sphereCentre[sphereIndex];
    for (unsigned short i=0; i<3; i++) {
        vMax[i] += spatialSize;
    }
    return {vMin, vMax};
}
/*
std::array<std::array<double, 2>, 3> Polysphere::getMinMaxVoxelCoordinates(size_t sphereIndex, const Vector<3> &position, const Orientation<3> &orientation, double spatialSize, const Orientation<3> &angularSize) const{
    std::array<double, 2>   sinAlpha = Polysphere::minmaxSin(orientation[0], angularSize[0]),
                            cosAlpha = Polysphere::minmaxCos(orientation[0], angularSize[0]),
                            sinBeta  = Polysphere::minmaxSin(orientation[1], angularSize[1]),
                            cosBeta  = Polysphere::minmaxCos(orientation[1], angularSize[1]),
                            sinGamma = Polysphere::minmaxSin(orientation[2], angularSize[2]),
                            cosGamma = Polysphere::minmaxCos(orientation[2], angularSize[2]);

    double xmin = position[0] +
        (
            Polysphere::sphereCentre[sphereIndex][0] * (cosAlpha[0] * cosGamma[0] - sinAlpha[1] * cosBeta[1] * sinGamma[1]) +
            Polysphere::sphereCentre[sphereIndex][1] * (-cosAlpha[1] * sinGamma[1] - sinAlpha[1] * cosBeta[1] * cosGamma[1]) +
            Polysphere::sphereCentre[sphereIndex][2] * (sinAlpha[0] * sinBeta[0])
        );
    double xmax = position[0] + spatialSize +
        (
            Polysphere::sphereCentre[sphereIndex][0] * (cosAlpha[1] * cosGamma[1] - sinAlpha[0] * cosBeta[0] * sinGamma[0]) +
            Polysphere::sphereCentre[sphereIndex][1] * (-cosAlpha[0] * sinGamma[0] - sinAlpha[0] * cosBeta[0] * cosGamma[0]) +
            Polysphere::sphereCentre[sphereIndex][2] * (sinAlpha[1] * sinBeta[1])
        );

    double ymin = position[1] +
        (
            Polysphere::sphereCentre[sphereIndex][0] * (sinAlpha[0] * cosGamma[0] + cosAlpha[0] * cosBeta[0] * sinGamma[0]) +
            Polysphere::sphereCentre[sphereIndex][1] * (- sinAlpha[1] * sinGamma[1] + cosAlpha[0] * cosBeta[0] * cosGamma[0]) +
            Polysphere::sphereCentre[sphereIndex][2] * (- cosAlpha[1] * sinBeta[1])
        );

    double ymax = position[1] + spatialSize +
        (
            Polysphere::sphereCentre[sphereIndex][0] * (sinAlpha[1] * cosGamma[1] + cosAlpha[1] * cosBeta[1] * sinGamma[1]) +
            Polysphere::sphereCentre[sphereIndex][1] * (- sinAlpha[0] * sinGamma[0] + cosAlpha[1] * cosBeta[1] * cosGamma[1]) +
            Polysphere::sphereCentre[sphereIndex][2] * (- cosAlpha[0] * sinBeta[0])
        );

    double zmin = position[2] +
        (
            Polysphere::sphereCentre[sphereIndex][0] * (sinBeta[0] * sinGamma[0]) +
            Polysphere::sphereCentre[sphereIndex][1] * (sinBeta[0] * cosGamma[0]) +
            Polysphere::sphereCentre[sphereIndex][2] * (cosBeta[0])
        );

    double zmax = position[2] + spatialSize +
        (
            Polysphere::sphereCentre[sphereIndex][0] * (sinBeta[1] * sinGamma[1]) +
            Polysphere::sphereCentre[sphereIndex][1] * (sinBeta[1] * cosGamma[1]) +
            Polysphere::sphereCentre[sphereIndex][2] * (cosBeta[1])
        );

    return std::array<std::array<double, 2>, 3>{std::array<double, 2>{xmin, xmax}, std::array<double, 2>{ymin, ymax}, std::array<double, 2>{zmin, zmax}};;
}
*/

bool Polysphere::voxelInside(BoundaryConditions<3> *bc, const RSAVector &voxelPosition, const RSAOrientation &voxelOrientation, double spatialSize, const RSAOrientation &angularSize) const{

    Orientation<3> angularVoxelSize = Shape<3, 3>::getAngularVoxelSize();
    if (
            (voxelOrientation[0] > angularVoxelSize[0]) ||
            (voxelOrientation[1] > angularVoxelSize[1]) ||
            (voxelOrientation[2] > angularVoxelSize[2])
    )
        return true;

    switch(this->voxelInsideEarlyRejection(bc, voxelPosition, voxelOrientation, spatialSize, angularSize)) {
        case TRUE:      return true;
        case FALSE:     return false;
        case UNKNOWN:   break;
    }

    // Version optimized for full angle - it was made mainly for OrientedFibrinogen, because it took ages to generate
    if (angularSize[0] >= 2*M_PI && angularSize[1] >= M_PI && angularSize[2] >= 2*M_PI)
        return this->fullAngleVoxelInside(bc, voxelPosition, spatialSize);

    Vector<3> translation = bc->getTranslation(voxelPosition, this->getPosition());
    Vector<3> thisPosition = this->getPosition() + translation;
    Orientation<3> thisOrientation = this->getOrientation();

    // loop over disks in virtual particle inside voxel
    for (size_t i = 0; i < Polysphere::sphereCentre.size(); i++){
        std::array<Vector<3>, 2> minMaxVoxelCoordinates = this->getMinMaxVoxelCoordinates(i, voxelPosition, voxelOrientation, spatialSize, angularSize);
        Vector<3> vMin = minMaxVoxelCoordinates[0];
        Vector<3> vMax = minMaxVoxelCoordinates[1];
        // loop over disks in this particle
        for (size_t j = 0; j < Polysphere::sphereCentre.size(); j++){
            // if a disk of virtual particle (placed anywhere in the voxel) intersects with disk of this particle, the voxel is fully inside the exclusion zone
            Vector<3> coordinates = Polysphere::getStaticSpherePosition(j, thisPosition, thisOrientation);

            double xmax2 = std::max( (vMin[0] - coordinates[0])*(vMin[0] - coordinates[0]) ,(vMax[0] - coordinates[0])*(vMax[0] - coordinates[0]));
            double ymax2 = std::max( (vMin[1] - coordinates[1])*(vMin[1] - coordinates[1]) ,(vMax[1] - coordinates[1])*(vMax[1] - coordinates[1]));
            double zmax2 = std::max( (vMin[2] - coordinates[2])*(vMin[2] - coordinates[2]) ,(vMax[2] - coordinates[2])*(vMax[2] - coordinates[2]));
            double r2 = (Polysphere::sphereR[i] + Polysphere::sphereR[j]) * (Polysphere::sphereR[i] + Polysphere::sphereR[j]);
            if (r2 > xmax2 + ymax2 + zmax2)
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
        Vector<3> spherePosition = this->getSpherePosition(i);
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
        Vector<3> spherePosition = this->getSpherePosition(i);
        out << "Sphere[{" << spherePosition[0] << ", " << spherePosition[1] << ", " << spherePosition[2] << "}, " << Polysphere::sphereR[i] << "]";
        if (i < max - 1)
            out << ", ";
    }
    out << "}";

    return out.str();
}