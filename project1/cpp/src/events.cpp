#include "../include/events.h"


// extern const double endTime;

Events events;

void run_events()
{
    double time = 0., dt = 0.001;
    int istep = 0;

    while (time <= endTime)
    {
        if (istep > 0) set_timestep(uf, dt, grid.delta, mu);

        events.runEvents(istep, time, dt);

        time += dt;
#if 0
        if (time > endTime && time - dt < endTime)
        {
            dt = endTime - (time - dt); // gets overridden by set_timestep. need to fix
            time = endTime;
        }
#endif
        istep++;
    }
}

