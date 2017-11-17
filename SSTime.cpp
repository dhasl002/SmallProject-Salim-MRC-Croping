
#include "SSTime.h"



void SSTime::startTimer()
{

    start =  std::chrono::system_clock::now();

}


double SSTime::getTimeElapsedSinceStart()
{
    std::chrono::time_point<std::chrono::system_clock> end =  std::chrono::system_clock::now();

    std::chrono::duration<double> elapsedSeconds = end-start;
    return elapsedSeconds.count();

}
