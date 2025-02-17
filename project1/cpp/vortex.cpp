#include "include/2D-grid.h"
#include "include/tracer.h"


const double width = 1;
const int level = 7;
const double endTime = 10;

#ifdef MU
const double mu = MU;
#else
const double mu = 0;
#endif


static void output_data (int istep, double time, double dt)
{
    if (istep == 0) {
        std::cout << "# of vector fields = " << VectorField::allVectorFields.size() << "\n";
        std::cout << "# of face vector fields = " << FaceVectorField::allFaceVectorFields.size() << "\n";
        std::cout << "# of scalar fields = " << ScalarField::allScalarFields.size() << "\n";
        std::cout << "# of events = " << events.allEvents.size() << "\n";
        std::cout << "\ni \t time \t dt\n";
    }
    
    std::cout << istep << " " << time << " " << dt << "\n";

    // output tracer and u field every 100 time steps
    if (istep % 1000 == 0) {
        std::string log1 = "data/" + std::to_string(istep) + "-tracer-" + std::to_string(time);
        std::ofstream output(log1);
        for (int i = ghost; i < grid.col - ghost; ++i) {
            for (int j = ghost; j < grid.row - ghost; ++j) {
                output << grid.x(i,j) << " " << grid.y(i,j) << " " << 
                       u.x(i,j) << " " << u.y(i,j) << " " << tracer(i,j) << "\n";
            }
            output << "\n";
        }
        output.close();
    }
}


SetEvent event(output_data);


int main() 
{
    init_grid (width, width, level);

    run_events();

    std::cout << "events done\nprinting solution to log\n";
    
    std::ofstream log;
    log.open("log");
    for (int i = ghost; i < grid.col - ghost; ++i) {
        for (int j = ghost; j < grid.row - ghost; ++j) {
            log << grid.x(i,j) << " " << grid.y(i,j) << " " << 
                   u.x(i,j) << " " << u.y(i,j) << " " << tracer(i,j) << "\n";
        }
        log << "\n";
    }
    log.close();
    
    return 0;
}
