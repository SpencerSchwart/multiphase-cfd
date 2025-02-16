#include "../../include/events.h"


// extern const double endTime;

Events events;


void run_events()
{
    double time = 0., dt = 0.1;
    int istep = 0;

    while (time < 1.)
    {
        events.runEvents(istep, time, dt);

        time += dt;
        istep++;
    }
}

