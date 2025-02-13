#ifndef 2D_GRID_H
#define 2D_GRID_H


#include <cstdlib>
#include <cstdio>
#include <functional>
#include <cmath>
#include <vector>
#include "timstep.h"
#include "math.h"


const int dimension = 2;
const int ghost = 1;    // layers of ghost cells around a domain boundary
double delta = 1.;


#define FOREACH_CELL() \
    for (int i = 0; i < grid.col; ++i) \
        for (int j = 0; j < grid.row; ++j)


#define FOREACH() \
    for (int i = ghost; i < grid.col - 2*ghost; ++i) \
        for (int j = ghost; j < grid.row - 2*ghost; ++j)


enum Boundaries {
    LEFT,           // 0
    RIGHT,          // 1
    TOP,            // 2
    BOTTOM,         // 3
    NUM_BOUNDARIES  // 4
};


struct Grid 
{
    double width, height;
    int level;

    double delta;
    int row, col;

    std::vector<double> x;
    std::vector<double> y;

    Grid() = default;

    Grid(double widthf, double heightf, int levelf)
        : heightf(height), width(width), level(levelf)
    {
        col(pow(2, level) + 2*ghost); row(pow(2, level) + 2*ghost);
        x((col * row), 0.0);
        y((col * row), 0.0);

        delta(max(width, height) / ((double) pow(2, levelf)));
    }
};


Grid grid;


struct Array2D
{
    int nx = 0, ny = 0;
    std::vector<double> data;
    
    Array2D() : 
        nx(grid.row), ny(grid.col), data(grid.width * grid.height, 0.0)
    {
    }

    Array2D(int rows, int cols, double initVal = 0.0)
        : nx(rows), ny(cols), data(rows * cols, initVal)
    {
    }

    void resize(int rows, int cols, double newVal = 0.0)
    {
        nx = rows; ny = cols;
        data.resize(nx * ny, newVal);
    }

    double& operator()(int i, int j)
    {
        return data[i * ny + j];
    }

    double operator()(int i, int j) const
    {
        return data[i * ny + j];
    }
};


struct Boundary
{
    std::array<std::function<double()>, NUM_BOUNDARIES> x;
    std::array<std::function<double()>, NUM_BOUNDARIES> y;

    bool dirichlet_x = false, dirichlet_y = false;
};


struct VectorField
{
    int nx = 0, ny = 0;
    Array2D x, y;

    VectorField() :
        nx(grid.row), ny(grid.col)
    {
        x(grid.row, grid.col, 0.0);
        y(grid.row, grid.col, 0.0);
    }

    VectorField(int rows, int cols, double initVal = 0.0)
        : x(rows, cols, initVal), y(rows, cols, initVal)
    {
    }

    void resize(int rows, int cols, int newVal = 0.0)
    {
        x.resize(rows, cols, newVal); y.resize(rows, cols, newVal);
    }

    Boundary boundary;
};


struct ScalarField
{
    int nx = 0, ny = 0;
    std::vector<double> data;

    ScalarField() :
        nx(grid.col), ny(grid.row), data(grid.col * grid.row, 0.0)
    {
    }

    ScalarField(int rows, int cols, double initVal = 0.0)
        nx(grid.col), ny(grid.row), data(grid.col * grid.row, initVal)
    {
    }

    void resize(int rows, int cols, double newVal = 0.0)
    {
        nx = rows; ny = cols;
        data.resize(nx * ny, newVal);
    }

    double& operator()(int i, int j)
    {
        return data[i * ny + j];
    }

    double operator()(int i, int j) const
    {
        return data[i * ny + j];
    }

    Boundary boundary;
};


void init_grid (double width, double height, int level);


bool is_ghost_cell (int i, int j);


void update_boundary_impl (VectorField& vf);


void update_boundary_impl (ScalarField& vf);


template <template... Fields>
void update_boundary (Fields&... fields) {
    (update_boundary_impl(fields), ...);
}


#endif // 2D_GRID_H
