//
// Created by PKua on 19.05.18.
//

#ifndef RSA3D_INFORMINGLOOPER_H
#define RSA3D_INFORMINGLOOPER_H


#include <iostream>
#include <string>

class InfoLooper {
private:
    std::ostream &out;
    std::size_t max;
    std::size_t infoStep;
    std::size_t i{};
    std::string text;

public:
    InfoLooper(size_t max, size_t infoStep, std::string text)
            : out(std::cout), max(max), infoStep(infoStep), text(std::move(text)) {}

    InfoLooper(size_t max, size_t infoStep, std::string text, std::ostream &out)
            : out(out), max(max), infoStep(infoStep), text(std::move(text)) {}

    bool step() {
        i++;
        if (i % infoStep == 0)
            out << ">> " << i << " " << text << std::endl;
        return i <= max;
    }

    size_t getIteration() const { return i; }
};


#endif //RSA3D_INFORMINGLOOPER_H
