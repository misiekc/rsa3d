//
// Created by PKua on 28.12.17.
//

#include "Timer.h"

using namespace std::chrono;

void Timer::start() {
    startTime = std::chrono::system_clock::now();
}

void Timer::stop() {
    endTime = std::chrono::system_clock::now();
}
