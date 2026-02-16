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
    const RND rnd(12345);
    size_t inside = 0;

    for (size_t i=0; i<mcTrials; i++){
        double x = (2.0 * rnd.nextValue() - 1.0) * maxSize;
        double y = (2.0 * rnd.nextValue() - 1.0) * maxSize;
        double z = (2.0 * rnd.nextValue() - 1.0) * maxSize;
        for(size_t index = 0; index < Polysphere::sphereR.size(); index++){
            Vector<3> coordinates = Polysphere::sphereCentre[i];
            double dx = coordinates[0] - x;
            double dy = coordinates[1] - y;
            double dz = coordinates[2] - z;
            if (dx*dx + dy*dy +dz*dz < Polysphere::sphereR[index]*Polysphere::sphereR[index]){
                inside++;
                break;
            }
        }
    }
    return (static_cast<double>(inside)/static_cast<double>(mcTrials))*8.0*maxSize*maxSize*maxSize;
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

std::array<std::array<double, 2>, 2> Polysphere::minmaxSinCos(double theta, double dt, double sinTheta, double cosTheta, double sinDt, double cosDt){

    theta = std::fmod(theta, 2*M_PI);
    if (theta<0)
        theta += 2*M_PI;

    double s1 = sinTheta;
    double s2 = sinTheta*cosDt + cosTheta*sinDt;

    double c1 = cosTheta;
    double c2 = cosTheta*cosDt - sinTheta*sinDt;

    std::array<double, 2> sr{
            std::min(s1, s2),
            std::max(s1,s2)
    };
    std::array<double, 2> cr{
            std::min(c1, c2),
            std::max(c1,c2)
    };

    if ( (theta < 0 && theta+dt > 0 ) || (theta < TWO_PI && theta+dt > TWO_PI) ){
        cr[1] = 1.0;
    }else if (theta < M_PI && theta+dt > M_PI){
        cr[0] = -1.0;
    }

    if (theta < HALF_PI && theta+dt > HALF_PI ){
        sr[1] = 1.0;
    }else if (theta < ONEHALF_PI && theta+dt > ONEHALF_PI){
        sr[0] = -1.0;
    }

    return {sr, cr};
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

std::array<Matrix<3,3>, 2> Polysphere::getMinMaxMatrices(const Orientation<3> &voxelOrientation, const Orientation<3> &angularSize) const {

    Orientation<3> translatedVoxelOrientation = Polysphere::fromVoxel(voxelOrientation);
    double angularSize1 = std::asin(voxelOrientation[1]+angularSize[1]-1.0) - std::asin(voxelOrientation[1]-1.0);

    Orientation<3> sinTheta{ std::sin(translatedVoxelOrientation[0]), std::sin(translatedVoxelOrientation[1]), std::sin(translatedVoxelOrientation[2]) };
    Orientation<3> cosTheta{ std::cos(translatedVoxelOrientation[0]), std::cos(translatedVoxelOrientation[1]), std::cos(translatedVoxelOrientation[2]) };
    Orientation<3> sinDt{ std::sin(angularSize[0]), std::sin(angularSize1), std::sin(angularSize[2])};
    Orientation<3> cosDt{ std::cos(angularSize[0]), std::cos(angularSize1), std::cos(angularSize[2])};

    std::array<std::array<double, 2>, 2> minmaxSinCos = Polysphere::minmaxSinCos(translatedVoxelOrientation[0], angularSize[0], sinTheta[0], cosTheta[0], sinDt[0], cosDt[0]);
    std::array<double, 2> sin_ax = minmaxSinCos[0];
    std::array<double, 2> cos_ax = minmaxSinCos[1];
    minmaxSinCos = Polysphere::minmaxSinCos(translatedVoxelOrientation[1], angularSize1, sinTheta[1], cosTheta[1], sinDt[1], cosDt[1]);
    std::array<double, 2> sin_ay = minmaxSinCos[0];
    std::array<double, 2> cos_ay = minmaxSinCos[1];
    minmaxSinCos = Polysphere::minmaxSinCos(translatedVoxelOrientation[2], angularSize[2], sinTheta[2], cosTheta[2], sinDt[2], cosDt[2]);
    std::array<double, 2> sin_az = minmaxSinCos[0];
    std::array<double, 2> cos_az = minmaxSinCos[1];

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
    return{minMatrix, maxMatrix};
}

std::array<Vector<3>, 2> Polysphere::getMinMaxVoxelCoordinates(size_t sphereIndex, const Vector<3> &voxelPosition, double spatialSize, const std::array<Matrix<3,3>, 2> &minmaxMatrices) const {
    Vector<3> vMin = voxelPosition + minmaxMatrices[0]*Polysphere::sphereCentre[sphereIndex];
    Vector<3> vMax = voxelPosition + minmaxMatrices[1]*Polysphere::sphereCentre[sphereIndex];
    for (unsigned short i=0; i<3; i++) {
        vMax[i] += spatialSize;
    }
    return {vMin, vMax};
}


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

    std::array<Matrix<3,3>, 2> minmaxMatrices = this->getMinMaxMatrices(voxelOrientation, angularSize);

    // loop over disks in virtual particle inside voxel
    for (size_t i = 0; i < Polysphere::sphereCentre.size(); i++){
        std::array<Vector<3>, 2> minMaxVoxelCoordinates = this->getMinMaxVoxelCoordinates(i, voxelPosition, spatialSize, minmaxMatrices);
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