//
// Created by PKua on 17.06.18.
//

#ifndef RSA3D_PARALELLINFOLOOPER_H
#define RSA3D_PARALELLINFOLOOPER_H


#include <iostream>
#include "../../rsa3d/Utils.h"


class ParallelInfoLooper {
private:
    std::ostream &out;
    std::size_t infoStep;
    std::size_t toReport;
    std::size_t counter{};
    std::string text;

public:
    explicit ParallelInfoLooper(size_t infoStep)
            : out(std::cout), infoStep(infoStep), toReport(infoStep), text(" pairs tested...") { }

    ParallelInfoLooper(size_t infoStep, std::string text)
            : out(std::cout), infoStep(infoStep), toReport(infoStep), text(std::move(text)) { }

    ParallelInfoLooper(size_t infoStep, std::string text, std::ostream &out)
            : out(out), infoStep(infoStep), toReport(infoStep), text(std::move(text)) { }

    void step() {
        _OMP_ATOMIC
        counter++;
        if (counter > toReport) {
            _OMP_CRITICAL(stdout)
            {
                if (counter > toReport) {
                    out << ">> " << toReport << text << std::endl;
                    toReport += infoStep;
                }
            }
        }
    }
};


#endif //RSA3D_PARALELLINFOLOOPER_H
