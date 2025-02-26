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

extern VectorField u;
void run_events()
{
    double time = 0., dt = 0.1;
    int istep = 0;

    while (time <= endTime)
    {
        //if (istep > 0) set_timestep(uf, dt, grid.delta, mu);
        std::cout << istep << " " << time << " " << dt << " " 
                  << Event::allEvents.size() << "\n";
        events.runEvents(istep, time, dt);
        stability (istep, time, dt);
        //std::cout << "top =" << u.boundary.x[TOP]() << "\n";
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

