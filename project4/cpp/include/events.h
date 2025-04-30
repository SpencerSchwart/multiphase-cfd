#ifndef EVENTS_H
#define EVENTS_H

#include <vector>
#include <functional>
#include "2D-grid.h"
#include "timestep.h"

extern FaceVectorField uf;
extern FaceVectorField mu;
extern double tend;

extern bool stop;

struct Event
{
    using Function = std::function<void(int,double,double&)>;

    Function func;
    std::string name;

    Event() = default;

    Event(const Function &userFunc, const std::string &userName, int index = -1)
    {
        func = userFunc;
        name = userName;
        if (index == -1)
            allEvents.push_back(this);
        else 
        {
            assert(index >= 0 
                && static_cast<size_t>(index) <= allEvents.size()
                && "ERROR: starting index for event is out of range!");
                
            allEvents.insert(allEvents.begin() + index, this);
        }
    }

    static std::vector<Event*> allEvents;
};

struct EventsManager
{   
    void runEvents(int i, double t, double dt);

    void printEvents();
};

extern EventsManager events;


void run_events();


#endif // EVENTS_H
