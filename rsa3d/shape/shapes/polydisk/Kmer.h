//
// Created by pkua on 15.11.2019.
//

#ifndef RSA3D_KMER_H
#define RSA3D_KMER_H


#include "Polydisk.h"

class Kmer : public Polydisk {
public:
    static void initClass(const std::string &attr);

    static double calculateArea(int numberOfDisks, double distance);
    static std::string preparePolydiskAttr(int numberOfDisks, double distance);
};


#endif //RSA3D_KMER_H
