#ifndef GRID_2D_H
#define GRID_2D_H

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <vector>
#include "math.h"

#define DCOL grid.col  // defualt # of columns
#define DROW grid.row  // defualt # of rows

const int dimension = 2;
const int ghost = 1;    // layers of ghost cells around a domain boundary
extern double delta;


#define FOREACH_CELL() \
    for (int i = 0; i < grid.col; ++i) \
        for (int j = 0; j < grid.row; ++j)


#define FOREACH() \
    for (int i = ghost; i < grid.col - ghost; ++i) \
        for (int j = ghost; j < grid.row - ghost; ++j)


#define FOREACH_FACE_X() \
    for (int i = ghost; i < grid.col - (ghost == 2); ++i) \
        for (int j = ghost; j < grid.row - ghost; ++j)

#define FOREACH_FACE_Y() \
    for (int i = ghost; i < grid.col - ghost; ++i) \
        for (int j = ghost; j < grid.row - (ghost == 2); ++j)


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
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j < ny);
        return data[i * ny + j];
    }

    double operator()(int i, int j) const
    {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j < ny);
        return data[i * ny + j];
    }
};


struct Grid 
{
    double width, height;
    int level;

    double delta;
    int row, col;

    Array2D x, y;

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

    std::array<bool, NUM_BOUNDARIES> dirichlet_x;
    std::array<bool, NUM_BOUNDARIES> dirichlet_y;
    
    Boundary()
    {
        for (int k = LEFT; k < NUM_BOUNDARIES; ++k)
        {
            dirichlet_x[k] = false;
            dirichlet_y[k] = false;
        }
    }
};


struct VectorField
{
    int nx = 0, ny = 0;
    Array2D x, y;

    VectorField() :
        nx(DROW), ny(DCOL)
    {
        x.resize(DROW, DCOL, 0.0);
        y.resize(DROW, DCOL, 0.0);

        allVectorFields.push_back(this);
    }

    VectorField(int rows, int cols, double initVal = 0.0)
        : x(rows, cols, initVal), y(rows, cols, initVal)
    {
        // allVectorFields.push_back(this);
    }

    void resize(int rows, int cols, int newVal = 0.0)
    {
        nx = cols; ny = rows;
        x.resize(rows, cols, newVal); y.resize(rows, cols, newVal);
    }

    Boundary boundary;

    static std::vector<VectorField*> allVectorFields;

};


struct FaceVectorField
{
    int nx = 0, ny = 0;
    Array2D x, y;

    FaceVectorField() :
        nx(DROW), ny(DCOL)
    {
        x.resize(DROW, DCOL, 0.0);
        y.resize(DROW, DCOL, 0.0);

        allFaceVectorFields.push_back(this);
    }

    FaceVectorField(int rows, int cols, double initVal = 0.0)
        : x(rows, cols, initVal), y(rows, cols, initVal)
    {
        allFaceVectorFields.push_back(this);
    }

    void resize(int rows, int cols, int newVal = 0.0)
    {
        nx = cols; ny = rows;
        x.resize(rows, cols, newVal); y.resize(rows, cols, newVal);
    }

    Boundary boundary;

    static std::vector<FaceVectorField*> allFaceVectorFields;
};


struct ScalarField
{
    int nx = 0, ny = 0;
    std::vector<double> data;

    ScalarField() :
        nx(DCOL), ny(DROW), data(DCOL * DROW, 0.0)
    {
        allScalarFields.push_back(this);
        for (int k = LEFT; k < NUM_BOUNDARIES; ++k)
        {
            dirichlet[k] = false;
        }
    }

    ScalarField(int rows, int cols, double initVal = 0.0)
        : nx(grid.col), ny(grid.row), data(grid.col * grid.row, initVal)
    {
        allScalarFields.push_back(this);
        for (int k = LEFT; k < NUM_BOUNDARIES; ++k)
        {
            dirichlet[k] = false;
        }
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

    static std::vector<ScalarField*> allScalarFields;

    std::array<std::function<double()>, NUM_BOUNDARIES> boundary;

    std::array<bool, NUM_BOUNDARIES> dirichlet;
};


void init_grid (double width, double height, int level);


bool is_ghost_cell (int i, int j);


void no_slip_boundary (VectorField& vf, Boundaries boundary, double val);


void no_slip_boundary (FaceVectorField& vf, Boundaries boundary, double val);


void outlet_boundary (ScalarField& sf, Boundaries boundary);


void update_boundary_impl (VectorField& vf);


void update_boundary_impl (ScalarField& vf);


void update_boundary_impl (FaceVectorField& vf);


template <typename... Fields>
void update_boundary (Fields&... fields) {
    (update_boundary_impl(fields), ...);
}


#include "events.h"
#endif // 2D_GRID_H
