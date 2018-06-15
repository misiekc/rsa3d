//
// Created by PKua on 31.03.18.
//

#include <memory>
#include <fstream>
#include <functional>
#include "AnisotropicShape2DExclusionDrawer.h"
#include "../rsa3d/Utils.h"
#include "../rsa3d/ShapeFactory.h"
#include "utility/MockBC.h"


namespace
{
    using namespace as2d_exdrawer;

    // Helper stuff declarations
    //------------------------------------------------------------------------------------------------------------------
    using overlap_funct = std::function<bool(const Vector<2> &)>;
    using shape_ptr = std::unique_ptr<Shape<2, 1>>;

    // Helper functions
    static Vector<2> bisect_ray(double initialRadius, double rayAngle, const overlap_funct &overlaps);
    static shape_ptr generate_shape(const Vector<2> &pos, double angle);
    static Polygon calculate_zone(double initialRadius, double resolution, const overlap_funct &overlapFunct);

    // Class containing basic options load from command line
    struct Context {
        std::string shape;
        std::string attr;
        double      shapeAngle{};
        double      angleFrom{};
        double      angleTo{};
        std::size_t resolution{};
        std::string outFilename{};

        void load(int argc, char **argv);
        shape_ptr generateMainShape();
    };


    // Helper stuff definitions
    //------------------------------------------------------------------------------------------------------------------
    void Context::load(int argc, char **argv) {
        if (argc != 9)
            die("Usage: rsa as2d_exdrawer <shape> <attr> <shape angle> <angle from> <angle to> <resolution> <output file>");
        this->shape = argv[2];
        this->attr = argv[3];
        this->shapeAngle = std::stod(argv[4]) * M_PI / 180.;
        this->angleFrom = std::stod(argv[5]) * M_PI / 180.;
        this->angleTo = std::stod(argv[6]) * M_PI / 180.;
        this->resolution = std::stoul(argv[7]);
        if (resolution == 0)
            die("Zero resolution");
        this->outFilename = argv[8];

        ShapeFactory::initShapeClass(this->shape, this->attr);
    }

    shape_ptr Context::generateMainShape() {
        return generate_shape(Vector<2>(), this->shapeAngle);
    }

    // Assuming overlaps is the characteristic function of a convex set, find its intersection with ary of angle
    // rayAngle by bisection starting from initialRadius
    Vector<2> bisect_ray(double initialRadius, double rayAngle, const overlap_funct &overlaps) {
        double radiusBeg = 0;
        double radiusEnd = initialRadius;
        Vector<2> normal{{std::cos(rayAngle), std::sin(rayAngle)}};

        // Initial radius should guarantee lack of overlap
        if (overlaps(normal * initialRadius))
            throw std::runtime_error("initialRadius overlap");

        // Bisect for the closest non-overlapping point until desired precision is reached
        double halfRadius = 0;
        while (radiusEnd - radiusBeg > 0.000001) {
            halfRadius = (radiusBeg + radiusEnd) / 2;
            if (overlaps(normal * halfRadius))
                radiusBeg = halfRadius;
            else
                radiusEnd = halfRadius;
        }
        return normal * halfRadius;
    }

    shape_ptr generate_shape(const Vector<2> &pos, double angle) {
        auto shape = dynamic_cast<Shape<2, 1>*>(ShapeFactory::createShape(nullptr));
        if (shape == nullptr)   die("Only Shape<2, 1> is supported");
        double arrayPos[] = {pos[0], pos[1]};
        shape->translate(arrayPos);
        shape->rotate({{angle}});
        return shape_ptr(shape);
    }

    // Return Polygon describing a boundary of convex set fo given characteristic function overlapFunct
    Polygon calculate_zone(double initialRadius, double resolution, const overlap_funct &overlapFunct) {
        Polygon result;

        double step = 2 * M_PI / resolution;
        double angle = 0;
        for (std::size_t i = 0; i < resolution; i++) {
            result.push_back(bisect_ray(initialRadius, angle, overlapFunct));
            angle += step;
        }

        return result;
    }

    // TODO duplicate code
    /* Uses pointInside method when dealing with ConvexShape or voxelInside with zero spatial size for normal Shape */
    bool point_inside(const Shape<2, 1> &shape, const Vector<2> &point, double angleFrom, double angleTo) {
        MockBC bc;
        double pointArr[2];
        point.copyToArray(pointArr);

        try {
            auto &convexShape = dynamic_cast<const ConvexShape<2, 1>&>(shape);
            return convexShape.pointInside(&bc, pointArr, {{angleFrom}}, angleTo - angleFrom);
        } catch (std::bad_cast&) {
            return shape.voxelInside(&bc, pointArr, {{angleFrom}}, 0, angleTo - angleFrom);
        }
    }
}


// Public interface definitions
//----------------------------------------------------------------------------------------------------------------------
namespace as2d_exdrawer
{
    int main(int argc, char **argv) {
        Context context;
        context.load(argc, argv);

        shape_ptr mainShape = context.generateMainShape();
        Polygon angleFromZone = zone_for_angle(*mainShape, context.angleFrom, context.resolution);
        Polygon angleToZone = zone_for_angle(*mainShape, context.angleTo, context.resolution);
        Polygon angleFromAndToZone = zone_for_two_angles(*mainShape, context.angleFrom,context.angleTo,
                                                         context.resolution);
        Polygon angleRangeZone = zone_for_angle_range(*mainShape, context.angleFrom, context.angleTo,
                                                      context.resolution);

        std::ofstream file(context.outFilename);
        if (!file)
            die("Cannot open file " + context.outFilename + " to write");
        print_notebook(*mainShape, angleFromZone, angleToZone, angleFromAndToZone, angleRangeZone, file);
        file.close();
        return EXIT_SUCCESS;
    }

    Polygon zone_for_angle(const Shape<2, 1> &shape, double angle, std::size_t resolution) {
        MockBC bc;
        overlap_funct overlapFunct = [&](const Vector<2> &pos) {
            auto secondShape = generate_shape(pos, angle);
            return shape.overlap(&bc, secondShape.get());
        };

        return calculate_zone(1.1 * Shape<2, 1>::getNeighbourListCellSize(), resolution, overlapFunct);
    }

    Polygon zone_for_two_angles(const Shape<2, 1> &shape, double angle1, double angle2, std::size_t resolution) {
        MockBC bc;
        overlap_funct overlapFunct = [&](const Vector<2> &pos) {
            auto secondShape1 = generate_shape(pos, angle1);
            auto secondShape2 = generate_shape(pos, angle2);
            return shape.overlap(&bc, secondShape1.get()) && shape.overlap(&bc, secondShape2.get());
        };

        return calculate_zone(1.1 * Shape<2, 1>::getNeighbourListCellSize(), resolution, overlapFunct);
    }

    Polygon zone_for_angle_range(const Shape<2, 1> &shape, double angleFrom, double angleTo,
                                 std::size_t resolution) {
        MockBC bc;
        overlap_funct overlapFunct = [&](const Vector<2> &pos) {
            double posArray[2];
            pos.copyToArray(posArray);
            return point_inside(shape, pos, angleFrom, angleTo);
        };

        return calculate_zone(1.1 * Shape<2, 1>::getNeighbourListCellSize(), resolution, overlapFunct);
    }

    void polygon_to_wolfram(const Polygon &polygon, std::ostream &stream) {
        stream << std::fixed;
        stream.precision(std::numeric_limits<double>::max_digits10);
        stream << "Polygon[{";
        for (auto vertex : Polygon(polygon.begin(), polygon.end() - 1)) {
            stream << vertex << "," << std::endl;
        }
        stream << polygon[polygon.size() - 1];
        stream << "}";
    }

    void print_notebook(const Shape<2, 1> &shape, const Polygon &fromZone, const Polygon &toZone,
                        const Polygon &fromAndToZone, const Polygon &rangeZone, std::ostream &stream) {
        stream << "(* Shape *)" << std::endl;
        stream << "shape = " << shape.toWolfram() << ";" << std::endl << std::endl;

        stream << "(* Angle from zone *)" << std::endl;
        stream << "fromZone = ";
        polygon_to_wolfram(fromZone, stream);
        stream << "];" << std::endl << std::endl;

        stream << "(* Angle to zone *)" << std::endl;
        stream << "toZone = ";
        polygon_to_wolfram(toZone, stream);
        stream << "];" << std::endl << std::endl;

        stream << "(* Angle from and to zone intersection *)" << std::endl;
        stream << "fromAndToZone = ";
        polygon_to_wolfram(fromAndToZone, stream);
        stream << "];" << std::endl << std::endl;

        stream << "(* Angle range zone *)" << std::endl;
        stream << "rangeZone = ";
        polygon_to_wolfram(rangeZone, stream);
        stream << "];" << std::endl << std::endl;

        stream << "Graphics[{Red, fromZone, Yellow, toZone, Purple, fromAndToZone, Orange, rangeZone, Black, shape}]";
        stream << std::endl << std::endl;
    }
}