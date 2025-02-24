#include "../include/events.h"

std::vector<Event*> Event::allEvents;

EventsManager events;


void EventsManager::runEvents(int i, double t, double dt)
{
    for (int _i = 0; _i < (int)Event::allEvents.size(); _i++)
    {
        Event::allEvents[_i]->func(i, t, dt);
#if 0
        std::cout << i << " " << t << " " << Event::allEvents[0]->name << "\n";
#endif
    }
}


void EventsManager::printEvents()
{
    for (int _i = 0; _i < (int)Event::allEvents.size(); _i++)
    {
        std::cout << _i << " " << Event::allEvents[_i]->name << "\n";
    }
}


void run_events()
{
    double time = 0., dt = 0.1;
    int istep = 0;

    while (time <= endTime)
    {
        //if (istep > 0) set_timestep(uf, dt, grid.delta, mu);
        std::cout << dt << " in event loop " << Event::allEvents.size() << "\n";
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

