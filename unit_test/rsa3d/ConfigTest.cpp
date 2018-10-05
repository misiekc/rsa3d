//
// Created by PKua on 27.09.18.
//

#include <catch.hpp>
#include <sstream>
#include "../../rsa3d/Config.h"

using namespace Catch::Matchers;

TEST_CASE("Config: basics and validation") {
    SECTION("empty") {
        std::istringstream str;

        auto config = Config::parse(str);

        REQUIRE(config.getKeys().empty());
    }

    SECTION("simple") {
        std::istringstream str("a=1\nb=2");

        auto config = Config::parse(str);

        REQUIRE(config.getKeys() == std::vector<std::string>{"a", "b"});
        REQUIRE(config.getString("a") == "1");
        REQUIRE(config.getString("b") == "2");
    }

    SECTION("custom delimiter") {
        std::istringstream str("a:1\nb:2");
        auto config = Config::parse(str, ':');

        REQUIRE(config.getKeys() == std::vector<std::string>{"a", "b"});
        REQUIRE(config.getString("a") == "1");
        REQUIRE(config.getString("b") == "2");
    }

    SECTION("# forbidden as delimiter") {
        std::istringstream str;

        REQUIRE_THROWS(Config::parse(str, '#'));
    }
}

TEST_CASE("Config: comments & whitespace") {
    SECTION("empty lines") {
        std::istringstream str("\na=1\n\nb=2\n");

        auto config = Config::parse(str);

        REQUIRE(config.getKeys() == std::vector<std::string>{"a", "b"});
        REQUIRE(config.getString("a") == "1");
        REQUIRE(config.getString("b") == "2");
    }

    SECTION("comment on a whole line") {
        std::istringstream str("#\n"
                               "#empty\n"
                               "b=2\n"
                               "#c=3commented field");

        auto config = Config::parse(str);

        REQUIRE(config.getKeys() == std::vector<std::string>{"b"});
        REQUIRE(config.getString("b") == "2");
    }

    SECTION("comments on a part of line") {
        std::istringstream str("a=1#comment appended\n"
                               "b=2\n");

        auto config = Config::parse(str);

        REQUIRE(config.getKeys() == std::vector<std::string>{"a", "b"});
        REQUIRE(config.getString("a") == "1");
        REQUIRE(config.getString("b") == "2");
    }

    SECTION("whitespace") {
        std::istringstream str("a= 1\nb =2\nc =  \nd\t=\t4");

        auto config = Config::parse(str);

        REQUIRE(config.getKeys() == std::vector<std::string>{"a", "b", "c", "d"});
        REQUIRE(config.getString("a") == "1");
        REQUIRE(config.getString("b") == "2");
        REQUIRE(config.getString("c").empty());
        REQUIRE(config.getString("d") == "4");
    }
}

TEST_CASE("Config: parsing errors") {
    SECTION("no delimiter") {
        std::istringstream str("a=1\nb2\n");

        REQUIRE_THROWS(Config::parse(str));
        // Somehow this breaks When(Method(...)).Return(...) in FakeIt (OMG WTF???)
        //REQUIRE_THROWS_WITH(Config::parse(str), Contains("line") && Contains("2"));
    }

    SECTION("field redefinition") {
        std::istringstream str("a=1\nb=2\na=3");

        REQUIRE_THROWS(Config::parse(str));
        //REQUIRE_THROWS_WITH(Config::parse(str), Contains("a") && Contains("line") && Contains("3"));
    }
}

TEST_CASE("Config: field access") {
    SECTION("non-existent") {
        std::istringstream str("a=1\nb=2\n");
        auto config = Config::parse(str);

        REQUIRE_THROWS(config.getString("c"));
        REQUIRE_THROWS(config.getInt("c"));
        REQUIRE_THROWS(config.getUnsignedLong("c"));
        REQUIRE_THROWS(config.getFloat("c"));
        REQUIRE_THROWS(config.getDouble("c"));
    }

    SECTION("correct") {
        std::istringstream str("string=foo\nint=-2\nuint=2\nfloat=-1.5e3\ndouble=-1.5e3");
        auto config = Config::parse(str);

        REQUIRE(config.getString("string") == "foo");
        REQUIRE(config.getInt("int") == -2);
        REQUIRE(config.getUnsignedLong("uint") == 2);
        REQUIRE(config.getFloat("float") == Approx(-1.5e3));
        REQUIRE(config.getDouble("double") == Approx(-1.5e3));
    }

    SECTION("incorrect") {
        std::istringstream str("minus=-1\ngarbage=foobar\n");
        auto config = Config::parse(str);

        REQUIRE_THROWS(config.getInt("garbage"));
        REQUIRE_THROWS(config.getUnsignedLong("minus"));
        REQUIRE_THROWS(config.getUnsignedLong("garbage"));
        REQUIRE_THROWS(config.getFloat("garbage"));
        REQUIRE_THROWS(config.getDouble("garbage"));
    }

    SECTION("hasField") {
        std::istringstream str("a=1\nb=2\n");
        auto config = Config::parse(str);

        REQUIRE(config.hasParam("a"));
        REQUIRE(config.hasParam("b"));
        REQUIRE(!config.hasParam("c"));
    }
}