#ifndef EVENTS_H
#define EVENTS_H

#include <vector>
#include <functional>



struct Events
{
    using Event = std::function<void(int,double,double)>;
    
    void addEvent (const Event& func)
    {
        allEvents.push_back(func);
    }

    void runEvents (int i, double t, double dt)
    {
        for (const auto& func : allEvents)
        {
            func(i, t, dt);
        }
    }

    std::vector<Event> allEvents;
};

extern Events events;

struct SetEvent
{
    using Event = std::function<void(int,double,double)>;

    SetEvent() = default;

    SetEvent(const Event& func)
    {
        events.addEvent(func);
    }
};

void run_events();


#endif // EVENTS_H
