//
// Created by ciesla on 15.11.2025.
//

#include "Kmer3D.h"

#include <assert.h>

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

std::array<Matrix<3,3>, 2> Kmer3D::getMinMaxMatrices(const Orientation<3> &translatedVoxelOrientation, const Orientation<3> &translatedAngularSize) const{

    Orientation<3> sinTheta{ std::sin(translatedVoxelOrientation[0]), std::sin(translatedVoxelOrientation[1]), std::sin(translatedVoxelOrientation[2]) };
    Orientation<3> cosTheta{ std::cos(translatedVoxelOrientation[0]), std::cos(translatedVoxelOrientation[1]), std::cos(translatedVoxelOrientation[2]) };
    Orientation<3> sinDt{ std::sin(translatedAngularSize[0]), std::sin(translatedAngularSize[1]), std::sin(translatedAngularSize[2])};
    Orientation<3> cosDt{ std::cos(translatedAngularSize[0]), std::cos(translatedAngularSize[1]), std::cos(translatedAngularSize[2])};


    std::array<std::array<double, 2>, 2> minmaxSinCos = Polysphere::minmaxSinCos(translatedVoxelOrientation[1], translatedAngularSize[1], sinTheta[1], cosTheta[1], sinDt[1], cosDt[1]);
    std::array<double, 2> sin_ay = minmaxSinCos[0];
    std::array<double, 2> cos_ay = minmaxSinCos[1];
    minmaxSinCos = Polysphere::minmaxSinCos(translatedVoxelOrientation[2], translatedAngularSize[2], sinTheta[2], cosTheta[2], sinDt[2], cosDt[2]);
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

// Normalize angle to [-pi, pi]
double normalize(double a)
{
    while (a <= -M_PI) a += 2*M_PI;
    while (a >   M_PI) a -= 2*M_PI;
    return a;
}

bool in_interval(double x, double a, double b)
{
    return x >= a && x <= b;
}

// Clamp angle into interval (used when optimum is outside)
double clamp(double x, double a, double b)
{
    return std::max(a, std::min(x, b));
}

// ---- cosine extrema on interval ----
double cos_min(double phi, double g0, double g1, double &arg)
{
    // candidate: cos = -1 at gamma = phi ± pi
    double g_star = normalize(phi + M_PI);

    if (in_interval(g_star, g0, g1)) {
        arg = g_star;
        return -1.0;
    }

    double c0 = std::cos(phi - g0);
    double c1 = std::cos(phi - g1);

    if (c0 < c1) {
        arg = g0;
        return c0;
    } else {
        arg = g1;
        return c1;
    }
}

double cos_max(double phi, double g0, double g1, double &arg)
{
    // candidate: cos = 1 at gamma = phi
    double g_star = normalize(phi);

    if (in_interval(g_star, g0, g1)) {
        arg = g_star;
        return 1.0;
    }

    double c0 = std::cos(phi - g0);
    double c1 = std::cos(phi - g1);

    if (c0 > c1) {
        arg = g0;
        return c0;
    } else {
        arg = g1;
        return c1;
    }
}

// ---- main solver ----
double minimize_expr(
    double A, const RSAVector &v, double b0, double b1, double g0, double g1, double &beta_out, double &gamma_out)
{
    double best = std::numeric_limits<double>::infinity();

    auto update = [&best, &beta_out, &gamma_out](double val, double b, double g){
        if (val < best) {
            best = val;
            beta_out = b;
            gamma_out = g;
        }
    };

    double R = std::sqrt(v[0]*v[0] + v[1]*v[1]);
    double phi = std::atan2(v[1], v[0]);
    double vnorm = v.norm();

    // ---- 1. interior stationary point ----
    if (vnorm > 0.0) {
        double beta_star = std::atan2(v[2], R);
        double gamma_star = phi;

        if (in_interval(beta_star, b0, b1) &&
            in_interval(gamma_star, g0, g1))
        {
            double val = -std::abs(A) * vnorm;
            update(val, beta_star, gamma_star);
        }
    }

    // ---- 2. beta = const edges ----
    for (double beta : {b0, b1})
    {
        double c = std::cos(beta);
        double s = std::sin(beta);

        double gamma_candidate;
        double cos_val;

        if (A * R * c >= 0)
            cos_val = cos_min(phi, g0, g1, gamma_candidate);
        else
            cos_val = cos_max(phi, g0, g1, gamma_candidate);

        double val = A * (R * c * cos_val - v[2] * s);
        update(val, beta, gamma_candidate);
    }

    // ---- 3. gamma = const edges ----
    for (double gamma : {g0, g1})
    {
        double M = v[0] * std::cos(gamma) + v[1] * std::sin(gamma);
        double amp = std::sqrt(M*M + v[2]*v[2]);

        if (amp == 0.0) continue;

        // ideal beta minimizing cos(...) = -1
        double beta_star = std::atan2(v[2], M) + M_PI;
        beta_star = normalize(beta_star);

        double beta_candidate;
        if (in_interval(beta_star, b0, b1))
            beta_candidate = beta_star;
        else
            beta_candidate = (std::abs(b0 - beta_star) < std::abs(b1 - beta_star)) ? b0 : b1;

        double val = A * (M * std::cos(beta_candidate) - v[2] * std::sin(beta_candidate));
        update(val, beta_candidate, gamma);
    }

    return best;
}

double maximize_expr(
    double A, const RSAVector &v, double b0, double b1, double g0, double g1, double &beta_out, double &gamma_out){
    // First find the minimizer of -f
    double beta_tmp, gamma_tmp;

    double min_neg = minimize_expr(-A, v, b0, b1, g0, g1, beta_tmp, gamma_tmp);

    // That point is the maximizer of original function
    beta_out = beta_tmp;
    gamma_out = gamma_tmp;

    // Evaluate the true maximum value
    double val = A * (std::cos(beta_out)*(v[0] * std::cos(gamma_out) + v[1] * std::sin(gamma_out))
                 - v[2] * std::sin(beta_out));

    return val;
}

double Kmer3D::voxelSphereMaxDistance(const RSAVector &voxelPosition, const RSAOrientation &translatedVoxelOrientation,
                                      double spatialSize, const RSAOrientation &translatedAngularSize,
                                      const RSAVector &virtualSphere, const RSAVector &thisSphere,
                                      RSAVector &vShapePosition, RSAOrientation &vShapeOrientation) {

    double dMax2 = -std::numeric_limits<double>::infinity();
    double virtualSphereLength = virtualSphere.norm();
    for (size_t ix=0; ix<=1; ix++) {
        for (size_t iy=0; iy<=1; iy++) {
            for (size_t iz=0; iz<=1; iz++) {
                RSAVector voxelVortex({voxelPosition[0]+ix*spatialSize, voxelPosition[1]+iy*spatialSize,
                                       voxelPosition[2]+iz*spatialSize});
                RSAVector voxelVector = voxelVortex - thisSphere;
                double voxelVectorLength = voxelVector.norm();
                double gamma_star, beta_star;

                double angularTerm = maximize_expr(virtualSphere[0], voxelVector,
                    translatedVoxelOrientation[1], translatedVoxelOrientation[1]+translatedAngularSize[1],
                    translatedVoxelOrientation[2], translatedVoxelOrientation[2]+translatedAngularSize[2],
                    beta_star,gamma_star);

                double dMaxVortex2 = voxelVectorLength*voxelVectorLength +
                            virtualSphereLength*virtualSphereLength +
                            2*angularTerm;
                
                if (dMaxVortex2 > dMax2) {
                    dMax2 = dMaxVortex2;
                    vShapeOrientation = {0, beta_star, gamma_star};
                    vShapePosition = voxelVortex;
                }
            }
        }
    }
    return std::sqrt(dMax2);
}

bool Kmer3D::isVoxelInside(BoundaryConditions<3> *bc, const RSAVector &voxelPosition, const RSAOrientation &translatedVoxelOrientation, double spatialSize, const RSAOrientation &translatedAngularSize) const{
    Vector<3> vShapePosition, vShapePositionMin;
    RSAOrientation vShapeOrientation, vShapeOrientationMin;
    int imax=-1, jmax=-1;
    Vector<3> thisPosition = this->getPosition();
    Orientation<3> thisOrientation = this->getOrientation();
    double dMax, dMinMax = std::numeric_limits<double>::infinity();
    Vector<3> voxelPositionLocal = voxelPosition + bc->getTranslation(thisPosition, voxelPosition);

    // loop over disks in virtual particle inside voxel
    for (size_t i = 0; i < Polysphere::sphereCentre.size(); i++){
        Vector<3> virtualSphere = Polysphere::sphereCentre[i];
        // loop over disks in this particle
        for (size_t j = 0; j < Polysphere::sphereCentre.size(); j++){
            Vector<3> thisSphere = Polysphere::getStaticSpherePosition(j, thisPosition, thisOrientation);
            dMax = Kmer3D::voxelSphereMaxDistance(voxelPositionLocal, translatedVoxelOrientation, spatialSize,
                                                         translatedAngularSize, virtualSphere, thisSphere,
                                                         vShapePosition, vShapeOrientation);
            if (Polysphere::sphereR[i] + Polysphere::sphereR[j] > dMax)
                return true;

            if (dMax < dMinMax) {
                dMinMax = dMax;
                vShapePositionMin = vShapePosition;
                vShapeOrientationMin = vShapeOrientation;
                imax = i;
                jmax = j;
            }
        }
    }
/*
    // before returning false we want to do an intersection test
    Shape *vs = new Kmer3D();
    vShapeOrientationMin[1] = std::sin(vShapeOrientationMin[1]) + 1.0;
    vs->rotate(vShapeOrientationMin);
    vShapeOrientationMin[1] = std::fmod(std::asin(vShapeOrientationMin[1] - 1.0)+2*M_PI, 2*M_PI);
    vs->translate(vShapePositionMin);
    if (this->overlap(bc, vs)) {

        this->overlap(bc, vs);

        int vIndex = 1, thisIndex = 1;

        std::cout << "wrong ";
        Vector<3> virtualSphere = Polysphere::sphereCentre[vIndex];
        double virtualSphereLength = virtualSphere.norm();
        Vector<3> thisSphere = Polysphere::getStaticSpherePosition(thisIndex, thisPosition, thisOrientation);
        RSAVector voxelVector = vs->getPosition() - thisSphere;
        double voxelVectorLength = voxelVector.norm();
        double gamma_star, beta_star;
        double angularTerm = maximize_expr(virtualSphere[0], voxelVector,
            translatedVoxelOrientation[1], translatedVoxelOrientation[1]+translatedAngularSize[1],
            translatedVoxelOrientation[2], translatedVoxelOrientation[2]+translatedAngularSize[2],
            beta_star,gamma_star);
        double dMaxVortex2 = voxelVectorLength*voxelVectorLength +
                    virtualSphereLength*virtualSphereLength +
                    2*angularTerm;

        this->overlap(bc, vs);
    }
    delete vs;
*/
    return false;
}

bool Kmer3D::voxelInside(BoundaryConditions<3> *bc, const Vector<3> &voxelPosition, const Orientation<3> &voxelOrientation, double spatialSize, const Orientation<3> &angularSize) const{
//    if (voxelOrientation[0] > 0.0)
    if (voxelOrientation[0] > 0.0 || voxelOrientation[1]>1.0)
        return true;
    else
        return Polysphere::voxelInside(bc, voxelPosition, voxelOrientation, spatialSize, angularSize);
}
