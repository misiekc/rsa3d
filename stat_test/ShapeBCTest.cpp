//
// Created by PKua on 12.08.18.
//

#include "ShapeBCTest.h"
#include "../rsa3d/utils/Utils.h"
#include "../rsa3d/shape/ShapeFactory.h"
#include "../rsa3d/boundary_conditions/FreeBC.h"
#include "../rsa3d/shape/ConvexShape.h"
#include "utils/UniformBallDistribution.h"
#include "utils/IndependentPairFactory.h"
#include "utils/InfoLooper.h"
#include "utils/TestExitCodes.h"

namespace {
    using ShapePair = RSAShapePairFactory::ShapePair;

    /* Conflicts tracker - counts them and stores example to display */
    struct Conflicts {
        std::size_t numOf{};
        ShapePair pair;
        ShapePair pairTranslated;
        RSAVector translation;

        void report(const ShapePair &pair, const ShapePair &pairTranslated, RSAVector translation);
        void print(std::ostream &out) const;
    };

    class Result {
    private:
        Conflicts ovConflicts;
        std::size_t overlapped{};
        Conflicts piConflicts;
        std::size_t pointInside{};

    public:
        void evaluatePair(const ShapePair &pair, const RSAVector &translation);
        void print(std::ostream &out) const;
        bool success() const { return ovConflicts.numOf == 0 && piConflicts.numOf == 0; }

        void evaluateOverlap(const ShapePair &pair, const ShapePair &pairTranslated, const RSAVector &translation);

        void evaluatePointInside(const ShapePair &pair, const ShapePair &pairTranslated, const RSAVector &translation);
    };

    /* BC returning programmed translation for a given pair. It expects positions of this shapes and returns a proper
     * translation - with a minus sign or not. If other points are tested, it throws. */
    class TranslatingBC : public RSABoundaryConditions {
    private:
        RSAVector translation;
        ShapePair testedPair;

        bool isFirstShape(const RSAVector &vector) const;
        bool isSecondShape(const RSAVector &vector) const;

    public:
        TranslatingBC(const RSAVector &translation, ShapePair testedPair) : translation(translation),
                                                                            testedPair(std::move(testedPair)) {}

        double distance2(const RSAVector &p1, const RSAVector &p2) const override;
        RSAVector getTranslation(const RSAVector &p1, const RSAVector &p2) const override;

    };

    bool TranslatingBC::isFirstShape(const RSAVector &vector) const {
        return (vector - testedPair.first()->getPosition()).norm2() < 1e-16;
    }

    bool TranslatingBC::isSecondShape(const RSAVector &vector) const {
        return (vector - testedPair.second()->getPosition()).norm2() < 1e-16;
    }

    double TranslatingBC::distance2(const RSAVector &p1, const RSAVector &p2) const {
        if (this->isFirstShape(p1) && this->isSecondShape(p2))          // translating p2 near p1
            return (p1 - p2 - translation).norm2();
        else if (this->isFirstShape(p2) && this->isSecondShape(p1))     // translating p1 near p2
            return (p1 - p2 + translation).norm2();
        else
            throw std::runtime_error("[ERROR] Unusual use of BoundaryConditions::distance2"
                                     " - unsupported by shape_bctest");
    }

    RSAVector TranslatingBC::getTranslation(const RSAVector &p1, const RSAVector &p2) const {
        if (this->isFirstShape(p1) && this->isSecondShape(p2))          // translating p2 near p1
            return translation;
        else if (this->isFirstShape(p2) && this->isSecondShape(p1))     // translating p1 near p2
            return -translation;
        else
            throw std::runtime_error("[ERROR] Unusual use of BoundaryConditions::getTranslation"
                                     " - unsupported by shape_bctest");
    }

    void Conflicts::report(const ShapePair &pair, const ShapePair &pairTranslated, RSAVector translation) {
        if (this->numOf == 0) {
            this->pair = pair;
            this->pairTranslated = pairTranslated;
            this->translation = translation;
        }
        this->numOf++;
    }

    void Conflicts::print(std::ostream &out) const {
        if (this->numOf == 0)   return;

        out << "[INFO] Evaluating the function above for this configuration: " << std::endl;
        this->pair.print(out);
        out << std::endl;
        out << "[INFO] using bc translation: " << this->translation << " yields different result" << std::endl;
        out << "[INFO] then for the same configuration with translation already applied:" << std::endl;
       this->pairTranslated.print(out);
    }

    void Result::print(std::ostream &out) const {
        out << "[INFO] " << overlapped << " pairs overlapped" << std::endl;

        if (ovConflicts.numOf == 0) {
            out << "[PASSED] Shape::overlap: No conflicts detected - BC are properly applied" << std::endl;
        } else {
            out << "[FAILED] Shape::overlap: " << ovConflicts.numOf << " conflicts detected" << std::endl;
            this->ovConflicts.print(out);
        }

        out << std::endl;
        out << "[INFO] " << pointInside << " shapes lie inside pointInside exclusion zone" << std::endl;

        // This part of a test will always be passed if Shape isn't ConvexShape
        if (piConflicts.numOf == 0) {
            out << "[PASSED] ConvexShape::pointInside: No conflicts detected - BC are properly applied" << std::endl;
        } else {
            out << "[FAILED] ConvexShape::pointInside: " << piConflicts.numOf << " conflicts detected" << std::endl;
            this->piConflicts.print(out);
        }
    }

    bool isConvexShape(const RSAShape *shape) {
        auto convexShape = dynamic_cast<const RSAConvexShape*>(shape);
        return convexShape != nullptr;
    }

    void Result::evaluatePair(const ShapePair &pair, const RSAVector &translation) {
        ShapePair pairTranslated{pair.first()->clone(), pair.second()->clone()};
        pairTranslated.second()->translate(translation);

        evaluateOverlap(pair, pairTranslated, translation);
        evaluatePointInside(pair, pairTranslated, translation);
    }

    void Result::evaluateOverlap(const ShapePair &pair, const ShapePair &pairTranslated, const RSAVector &translation) {
        TranslatingBC translatingBC{translation, pair};
        RSAFreeBC mockBC;

        bool bcOverlap = pair.first()->overlap(&translatingBC, pair.second());
        bool transOverlap = pairTranslated.first()->overlap(&mockBC, pairTranslated.second());

        if (bcOverlap)
            overlapped++;
        if (bcOverlap != transOverlap)
            ovConflicts.report(pair, pairTranslated, translation);
    }

    void Result::evaluatePointInside(const ShapePair &pair, const ShapePair &pairTranslated,
                                     const RSAVector &translation) {
        TranslatingBC translatingBC{translation, pair};
        RSAFreeBC mockBC;

        bool bcPointInside;
        bool transPointInside;

        if (isConvexShape(pair.first())) {
            auto &firstConvex = dynamic_cast<const RSAConvexShape &>(*pair.first());
            auto &secondConvex = dynamic_cast<const RSAConvexShape &>(*pair.second());
            auto &firstConvexTrans = dynamic_cast<const RSAConvexShape &>(*pairTranslated.first());
            auto &secondConvexTrans = dynamic_cast<const RSAConvexShape &>(*pairTranslated.second());

            bcPointInside = firstConvex.pointInside(&translatingBC, secondConvex.getPosition());
            transPointInside = firstConvexTrans.pointInside(&mockBC, secondConvexTrans.getPosition());
        } else {
            RSAOrientation zeroOrientation;
            zeroOrientation.fill(0);

            bcPointInside = pair.first()->voxelInside(&translatingBC, pair.second()->getPosition(), zeroOrientation, 0,
                    RSAShape::getAngularVoxelSize());
            transPointInside = pairTranslated.first()->voxelInside(&mockBC, pairTranslated.second()->getPosition(),
                    zeroOrientation, 0, RSAShape::getAngularVoxelSize());
        }

        if (bcPointInside)
            pointInside++;
        if (bcPointInside != transPointInside)
            piConflicts.report(pair, pairTranslated, translation);
    }

    Result perform_test(IndependentPairFactory &factory, unsigned long tries) {
        std::cout << "[INFO] Performing test for " << tries << " pairs using " << factory.getDescription() << std::endl;

        Result result;
        InfoLooper looper{tries, 10000, "pairs tested..."};
        Distribution &distribution = factory.getDistribution();
        RND rnd;
        while (looper.step())
            result.evaluatePair(factory.generate(), distribution.randomVector(&rnd));

        return result;
    }
}

namespace shape_bctest
{
    int main(int argc, char **argv)
    {
        if (argc < 6) {
            std::cerr << "Usage: ./rsa shape_bctest [particle] [attr] [ball_radius] [max_tries]" << std::endl;
            return TEST_ERROR;
        }

        double ballRadius = std::stod(argv[4]);
        unsigned long maxTries = std::stoul(argv[5]);
        if (ballRadius <= 0 || maxTries <= 0) {
            std::cerr << "Wrong input. Aborting." << std::endl;
            return TEST_ERROR;
        }

        ShapeFactory::initShapeClass(argv[2], argv[3]);
        RSAShape::setEarlyRejectionEnabled(false);
        UniformBallDistribution distribution{ballRadius};
        IndependentPairFactory factory{distribution};

        Result results = perform_test(factory, maxTries);
        std::cout << std::endl;
        results.print(std::cout);

        return results.success() ? TEST_SUCCESS : TEST_FAILURE;
    }
}