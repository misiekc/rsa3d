//
// Created by PKua on 03.07.18.
//

#ifndef RSA3D_ORDERCALCULABLE_H
#define RSA3D_ORDERCALCULABLE_H

#include <vector>

class OrderCalculable {
public:
    virtual ~OrderCalculable() = default;

    virtual std::vector<double> calculateOrder(const OrderCalculable *other) const = 0;
    virtual std::size_t getNumOfOrderParameters() const { return this->calculateOrder(this).size(); };
};


#endif //RSA3D_ORDERCALCULABLE_H
