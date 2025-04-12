#include "../include/events.h"

std::vector<Event*> Event::allEvents;

EventsManager events;

bool stop = false;


void EventsManager::runEvents(int i, double t, double dt)
{
    for (int _i = 0; _i < (int)Event::allEvents.size(); _i++)
        Event::allEvents[_i]->func(i, t, dt);
}


void EventsManager::printEvents()
{
    for (int _i = 0; _i < (int)Event::allEvents.size(); _i++)
        std::cout << _i << " " << Event::allEvents[_i]->name << "\n";
}


extern VectorField u;
void run_events()
{
    double time = 0., dt = dtmax;
    int istep = 0;

    while (time <= tend && !stop)
    {
        events.runEvents(istep, time, dt);
        stability (istep, time, dt);
        time += dt;

        if(time >= tend)
            stop = true;
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

