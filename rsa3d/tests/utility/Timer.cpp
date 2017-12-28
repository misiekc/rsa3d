//
// Created by PKua on 28.12.17.
//

#include "../CuboidSpeedTest.h"
#include <bits/unique_ptr.h>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include "Timer.h"

using namespace std::chrono;

void Timer::start() {
    startTime = std::chrono::_V2::system_clock::now();
}

void Timer::stop() {
    endTime = std::chrono::_V2::system_clock::now();
}
