//
// Created by PKua on 23.09.18.
//

//#####################//
// ENABLED FOR rsa.2.0 //
//#####################//

#if RSA_SPATIAL_DIMENSION == 2 && RSA_ANGULAR_DIMENSION == 0


#include "catch.hpp"
#include "fakeit.hpp"
#include "../../rsa3d/Packing.h"

using namespace fakeit;

TEST_CASE("Packing: ctor & dtor") {
    SECTION("default") {
        Packing packing;

        REQUIRE(packing.size() == 0);
    }

    SECTION("copy") {
        Mock<Shape<2, 0>> s1cloneMock;      Shape<2, 0> *s1clone = &s1cloneMock.get();
        Mock<Shape<2, 0>> s1mock;
        Method(s1mock,clone) = s1clone;
        Shape<2, 0> *s1 = &s1mock.get();

        Mock<Shape<2, 0>> s2cloneMock;      Shape<2, 0> *s2clone = &s2cloneMock.get();
        Mock<Shape<2, 0>> s2mock;
        Method(s2mock,clone) = s2clone;
        Shape<2, 0> *s2 = &s2mock.get();

        Packing packing;
        packing.addShape(s1);
        packing.addShape(s2);


        Packing packing2(packing);


        REQUIRE(packing.size() == 2);
        REQUIRE(packing2.size() == 2);
        REQUIRE(packing2[0] == s1clone);
        REQUIRE(packing2[1] == s2clone);
        Verify(Method(s1mock,clone)).Once();
        Verify(Method(s2mock,clone)).Once();
        VerifyNoOtherInvocations(s1mock, s2mock, s1cloneMock, s2cloneMock);
    }

    SECTION("assignment") {
        Mock<Shape<2, 0>> prevS1Mock;
        Fake(Dtor(prevS1Mock));
        Shape<2, 0> *prevS1 = &prevS1Mock.get();

        Mock<Shape<2, 0>> prevS2Mock;
        Fake(Dtor(prevS2Mock));
        Shape<2, 0> *prevS2 = &prevS2Mock.get();

        Mock<Shape<2, 0>> newSCloneMock;    Shape<2, 0> *newSClone = &newSCloneMock.get();
        Mock<Shape<2, 0>> newSMock;
        Method(newSMock,clone) = newSClone;
        Shape<2, 0> *newS = &newSMock.get();

        Packing packing;
        packing.addShape(newS);
        Packing packing2;
        packing2.addShape(prevS1);
        packing2.addShape(prevS2);


        packing2 = packing;


        REQUIRE(packing2.size() == 1);
        REQUIRE(packing2[0] == newSClone);
        Verify(Dtor(prevS1Mock));
        Verify(Dtor(prevS2Mock));
        Verify(Method(newSMock,clone)).Once();
        VerifyNoOtherInvocations(prevS1Mock, prevS2Mock, newSMock, newSCloneMock);
    }

    SECTION("move") {
        Mock<Shape<2, 0>> s1mock;   Shape<2, 0> *s1 = &s1mock.get();
        Mock<Shape<2, 0>> s2mock;   Shape<2, 0> *s2 = &s2mock.get();

        Packing packing;
        packing.addShape(s1);
        packing.addShape(s2);

        Packing packing2(std::move(packing));

        REQUIRE(packing.size() == 0);
        REQUIRE(packing2.size() == 2);
        REQUIRE(packing2[0] == s1);
        REQUIRE(packing2[1] == s2);
        VerifyNoOtherInvocations(s1mock, s2mock);
    }

    SECTION("move assignment") {
        Mock<Shape<2, 0>> prevS1Mock;
        Fake(Dtor(prevS1Mock));
        Shape<2, 0> *prevS1 = &prevS1Mock.get();

        Mock<Shape<2, 0>> prevS2Mock;
        Fake(Dtor(prevS2Mock));
        Shape<2, 0> *prevS2 = &prevS2Mock.get();

        Mock<Shape<2, 0>> newSMock; Shape<2, 0> *newS = &newSMock.get();

        Packing packing;
        packing.addShape(newS);
        Packing packing2;
        packing2.addShape(prevS1);
        packing2.addShape(prevS2);


        packing2 = std::move(packing);


        REQUIRE(packing2.size() == 1);
        REQUIRE(packing2[0] == newS);
        Verify(Dtor(prevS1Mock));
        Verify(Dtor(prevS2Mock));
        VerifyNoOtherInvocations(prevS1Mock, prevS2Mock, newSMock);
    }

    SECTION("destructor") {
        Mock<Shape<2, 0>> s1mock;
        Fake(Dtor(s1mock));
        Shape<2, 0> *s1 = &s1mock.get();
        Mock<Shape<2, 0>> s2mock;
        Fake(Dtor(s2mock));
        Shape<2, 0> *s2 = &s2mock.get();
        auto packing = new Packing();
        packing->addShape(s1);
        packing->addShape(s2);

        delete packing;

        Verify(Dtor(s1mock));
        Verify(Dtor(s2mock));
        VerifyNoOtherInvocations(s1mock, s2mock);
    }
}

TEST_CASE("Packing: add & remove") {
    SECTION("add") {
        Mock<Shape<2, 0>> s1mock;   Shape<2, 0> *s1 = &s1mock.get();
        Mock<Shape<2, 0>> s2mock;   Shape<2, 0> *s2 = &s2mock.get();
        Packing packing;

        packing.addShape(s1);
        packing.addShape(s2);

        REQUIRE(packing.size() == 2);
        REQUIRE(packing[0] == s1);
        REQUIRE(packing[1] == s2);
        VerifyNoOtherInvocations(s1mock, s2mock);
    }

    SECTION("remove") {
        Mock<Shape<2, 0>> s1mock;
        Fake(Dtor(s1mock));
        Shape<2, 0> *s1 = &s1mock.get();
        Mock<Shape<2, 0>> s2mock;   Shape<2, 0> *s2 = &s2mock.get();
        Packing packing;
        packing.addShape(s1);
        packing.addShape(s2);

        packing.removeShape(0);

        REQUIRE(packing.size() == 1);
        REQUIRE(packing[0] == s2);
        Verify(Dtor(s1mock));
        VerifyNoOtherInvocations(s1mock, s2mock);
    }

    SECTION("remove nonexistent") {
        Mock<Shape<2, 0>> sMock;   Shape<2, 0> *s = &sMock.get();
        Packing packing;
        packing.addShape(s);

        REQUIRE_THROWS(packing.removeShape(1));

        REQUIRE(packing.size() == 1);
        REQUIRE(packing[0] == s);
        VerifyNoOtherInvocations(sMock);
    }

    SECTION("clean") {
        Mock<Shape<2, 0>> s1mock;
        Fake(Dtor(s1mock));
        Shape<2, 0> *s1 = &s1mock.get();
        Mock<Shape<2, 0>> s2mock;
        Fake(Dtor(s2mock));
        Shape<2, 0> *s2 = &s2mock.get();
        Packing packing;
        packing.addShape(s1);
        packing.addShape(s2);

        packing.clear();

        REQUIRE(packing.size() == 0);
        Verify(Dtor(s1mock));
        Verify(Dtor(s2mock));
        VerifyNoOtherInvocations(s1mock, s2mock);
    }
}

TEST_CASE("Packing: access") {
    SECTION("non-existing") {
        Mock<Shape<2, 0>> sMock;   Shape<2, 0> *s = &sMock.get();
        Packing packing;
        packing.addShape(s);

        REQUIRE_THROWS(packing[1]);
        REQUIRE(packing.size() == 1);
        REQUIRE(packing[0] == s);
        VerifyNoOtherInvocations(sMock);
    }

    SECTION("back & front") {
        Mock<Shape<2, 0>> s1mock;   Shape<2, 0> *s1 = &s1mock.get();
        Mock<Shape<2, 0>> s2mock;   Shape<2, 0> *s2 = &s2mock.get();
        Packing packing;
        packing.addShape(s1);
        packing.addShape(s2);

        REQUIRE(packing.front() == s1);
        REQUIRE(packing.back() == s2);
        VerifyNoOtherInvocations(s1mock, s2mock);
    }

    SECTION("vector view") {
        Mock<Shape<2, 0>> s1mock;   Shape<2, 0> *s1 = &s1mock.get();
        Mock<Shape<2, 0>> s2mock;   Shape<2, 0> *s2 = &s2mock.get();
        Packing packing;
        packing.addShape(s1);
        packing.addShape(s2);

        auto vec = packing.getVector();

        REQUIRE(vec.size() == 2);
        REQUIRE(vec[0] == s1);
        REQUIRE(vec[1] == s2);
        VerifyNoOtherInvocations(s1mock, s2mock);
    }
}

TEST_CASE("Packing: PBC expand") {
    SECTION("top-right expand") {

        /*
         * s1 on the boundary, s2 not to be expanded
         *
         *    0   0.2  0.8   1
         * 0  +--------------+              +--------------+
         *    |           s1 |           s3 |           s1 |
         * 0.2|    +----+    |              |    +----+    |
         *    |    | s2 |    |   ----->     |    | s2 |    |
         * 0.8|    +----+    |              |    +----+    |
         *    |              |              |              |
         * 1  +--------------+              +--------------+
         *                               s5             s4
         */

        // Temporarily commented until issue with mocking resolved

        /*Mock<Shape<2, 0>> s4mock;
        Fake(Method(s4mock,translate));
        Shape<2, 0> *s4 = &s4mock.get();

        Mock<Shape<2, 0>> s5mock;
        Fake(Method(s5mock,translate));
        Shape<2, 0> *s5 = &s5mock.get();

        Mock<Shape<2, 0>> s3mock;
        Fake(Method(s3mock,translate));
        Method(s3mock,getPosition) = {{-0.1, 0.1}};
        Method(s3mock,clone) = s5;
        Shape<2, 0> *s3 = &s3mock.get();

        Mock<Shape<2, 0>> s1mock;
        Method(s1mock,getPosition) = {{0.9, 0.1}};
        When(Method(s1mock,clone)).Return(s3).Return(s4);
        Shape<2, 0> *s1 = &s1mock.get();

        Mock<Shape<2, 0>> s2mock;
        Method(s2mock,getPosition) = {{0.3, 0.3}};
        Shape<2, 0> *s2 = &s2mock.get();

        Packing packing;
        packing.addShape(s1);
        packing.addShape(s2);


        packing.expandOnPBC(1, 0.2);


        REQUIRE(packing.size() == 5);
        REQUIRE(packing.getVector() == std::vector<const Shape<2, 0>*>{s1, s2, s3, s4, s5});
        Verify(Method(s1mock,getPosition) + Method(s1mock,clone)).Twice();
        Verify(Method(s2mock,getPosition)).Twice();
        Verify(Method(s3mock,translate).Using({{-1, 0}}) + Method(s3mock,getPosition) + Method(s3mock,clone)).Once();
        Verify(Method(s4mock,translate).Using({{0, 1}})).Once();
        Verify(Method(s5mock,translate).Using({{0, 1}})).Once();
        VerifyNoOtherInvocations(s1mock, s2mock, s3mock, s4mock, s5mock);*/
    }

    SECTION("invalid arguments") {
        Packing packing;

        REQUIRE_THROWS(packing.expandOnPBC(0));
        REQUIRE_THROWS(packing.expandOnPBC(-1));
    }
}


#endif