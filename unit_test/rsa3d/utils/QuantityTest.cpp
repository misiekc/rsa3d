//
// Created by pkua on 06.06.19.
//

#include <catch.hpp>
#include <sstream>
#include "../../../rsa3d/utils/Quantity.h"

TEST_CASE("Quantity: rounding based on error") {
    SECTION("basic") {
        SECTION("decimal with error with fractional part") {
            std::ostringstream result1, result2;

            result1 << Quantity(34.344, 2.559, 3);          // log(x) = integer + less than 0.5
            result2 << Quantity(234.234544, 0.03425, 3);    // log(x) = integer + more than 0.5

            REQUIRE(result1.str() == "34.34\t2.56");
            REQUIRE(result2.str() == "234.2345\t0.0343");
        }

        SECTION("negative decimal") {
            std::ostringstream result;

            result << Quantity(-34.344, 2.559, 3);

            REQUIRE(result.str() == "-34.34\t2.56");
        }

        SECTION("decimal with error touching fractional part") {
            std::ostringstream result;

            result << Quantity(1345.34, 127.556, 3);

            REQUIRE(result.str() == "1345\t128");
        }

        SECTION("decimal with error not touching fractional part") {
            std::ostringstream result;

            result << Quantity(28345.34, 1277.556, 3);

            REQUIRE(result.str() == "28350\t1280");
        }

        /*SECTION("wrong 12.5 rounding") {
            std::ostringstream result;

            result << Quantity(2.343214, 0.125, 2);

            REQUIRE(result.str() == "2.34\t0.13");
        }*/
    }

    SECTION("tricky") {
        SECTION("error first significant digit changing after rounding") {
            std::ostringstream result;

            result << Quantity(0.294838, 0.009977, 2);

            REQUIRE(result.str() == "0.295\t0.010");
        }

        SECTION("decimal with larger error") {
            std::ostringstream result1, result2;

            result1 << Quantity(445.6, 2256.5, 3);
            result2 << Quantity(0.4563, 1.567, 3);

            REQUIRE(result1.str() == "450\t2260");
            REQUIRE(result2.str() == "0.46\t1.57");
        }

        SECTION("decimal with larger error and order changing after rounding") {
            std::ostringstream result;

            result << Quantity(0.999838, 1122.5, 2);

            REQUIRE(result.str() == "1\t1100");
        }

        SECTION("decimal with larger error before rounding") {
            std::ostringstream result;

            result << Quantity(0.999838, 12.65, 2);

            REQUIRE(result.str() == "1\t13");
        }

        SECTION("decimal with much larger error") {
            std::ostringstream result;

            result << Quantity(54.62, 256332.5, 3);

            REQUIRE(result.str() == "50\t256000");
        }

        SECTION("decimal with zeros in value") {
            std::ostringstream result;

            result << Quantity(45.6, 0.0002332, 3);

            REQUIRE(result.str() == "45.600000\t0.000233");
        }

        SECTION("decimal with zeros in error") {
            std::ostringstream result;

            result << Quantity(45.6, 0.0002, 3);

            REQUIRE(result.str() == "45.600000\t0.000200");
        }
    }
}