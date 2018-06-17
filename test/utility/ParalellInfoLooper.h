//
// Created by PKua on 17.06.18.
//

#ifndef RSA3D_PARALELLINFOLOOPER_H
#define RSA3D_PARALELLINFOLOOPER_H


#include <iostream>


#define STRINGIFY(x) #x

#ifdef _OPENMP
    #define _OMP_PARALLEL_FOR   _Pragma("omp parallel for")
    #define _OMP_ATOMIC         _Pragma("omp atomic")
    #define _OMP_CRITICAL(x)    _Pragma(STRINGIFY(omp critical(x)))
#else
#define _OMP_PARALLEL_FOR
    #define _OMP_ATOMIC
    #define _OMP_CRITICAL(x)
#endif


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
