#ifndef ____SSTime__
#define ____SSTime__

#include <stdio.h>
#include <iostream>
#include <chrono>
#include <ctime>


class SSTime {

private:
    std::chrono::time_point<std::chrono::system_clock> start;

public:
    SSTime () {

        start =  std::chrono::system_clock::now();
    };

    void startTimer();
    double getTimeElapsedSinceStart();

    ~SSTime () {
    };

};



#endif /* defined(____SSTime__) */
