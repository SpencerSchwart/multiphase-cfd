#ifndef GRID_2D_H
#define GRID_2D_H


#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <functional>
#include <cmath>
#include <vector>
//#include "timestep.h"
//#include "math.h"


const int dimension = 2;
const int ghost = 1;    // layers of ghost cells around a domain boundary


#define FOREACH_CELL() \
    for (int i = 0; i < grid.col; ++i) \
        for (int j = 0; j < grid.row; ++j)


#define FOREACH() \
    for (int i = ghost; i < grid.col - ghost; ++i) \
        for (int j = ghost; j < grid.row - ghost; ++j)


enum Boundaries {
    LEFT,           // 0
    RIGHT,          // 1
    TOP,            // 2
    BOTTOM,         // 3
    NUM_BOUNDARIES  // 4
};


struct Array2D
{
    int nx = 0, ny = 0;
    std::vector<double> data;
   
    Array2D() = default;

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


struct Grid 
{
    double width, height;
    int level;

    double delta;
    int row, col;

    Array2D x;
    Array2D y;

    Grid() = default;

    Grid(double widthf, double heightf, int levelf)
        : width(widthf), height(heightf), level(levelf)
    {
        int size = pow(2, level) + 2*ghost;
        col = size; row = size;
        x.resize((col * row), 0.0);
        y.resize((col * row), 0.0);

        delta = fmax(width, height) / ((double) pow(2, levelf));
    }

    void resize (double widthf, double heightf, int levelf)
    {
        width = widthf, height = heightf, level = levelf;
        int size = pow(2, level) + 2*ghost;
        col = size; row = size;
        x.resize(row, col, 0.0);
        y.resize(row, col, 0.0);

        delta = fmax(width, height) / ((double) pow(2, levelf));
    }
};


extern Grid grid;


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
        x.resize(grid.row, grid.col, 0.0);
        y.resize(grid.row, grid.col, 0.0);
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
        : nx(grid.col), ny(grid.row), data(grid.col * grid.row, initVal)
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


template <typename... Fields>
void update_boundary (Fields&... fields) {
    (update_boundary_impl(fields), ...);
}


#endif // 2D_GRID_H
