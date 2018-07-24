//
// Created by PKua on 15.07.18.
//

#ifndef RSA3D_SNUBDODECAHEDRON_H
#define RSA3D_SNUBDODECAHEDRON_H


#include "RegularSolid.h"
#include "../../OrderCalculable.h"

class SnubDodecahedron : public RegularSolid<SnubDodecahedron>, public OrderCalculable {
private:
    friend RegularSolid<SnubDodecahedron>;

    const static TriTriOverlapRS overlapStrategy;

    static void calculateStatic(const std::string &attr);

public:
    explicit SnubDodecahedron(const Matrix<3, 3> &orientation) : RegularSolid(orientation) {}

    OverlapStrategy<3, 0> *createStrategy(const std::string &name) const override;

    std::vector<double> calculateOrder(const OrderCalculable *other) const override;
};


#endif //RSA3D_SNUBDODECAHEDRON_H
