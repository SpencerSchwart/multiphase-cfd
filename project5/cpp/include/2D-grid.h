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

#define DCOL (grid.col)  // defualt # of columns
#define DROW (grid.row)  // defualt # of rows

#define NODATA (1e30)
#define NZERO  (1e-30)   // to avoid dividing by zero

const int dimension = 2;
const int ghost = 1;    // layers of ghost cells around a domain boundary
extern double& delta;

#define dv() (sq(delta))

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

#define FOREACH_VERTEX() \
    for (int i = ghost; i < grid.col - (ghost == 2); ++i) \
        for (int j = ghost; j < grid.row - (ghost == 2); ++j)


enum Boundaries 
{
    LEFT,           // 0
    RIGHT,          // 1
    TOP,            // 2
    BOTTOM,         // 3
    NUM_BOUNDARIES  // 4
};


struct Point
{
    int i, j;
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

    Array2D x, y;     // cell center coordinates
    Array2D xv, yv;   // vertex coordinates

    Grid() = default;

    Grid(double widthf, double heightf, int levelf)
        : width(widthf), height(heightf), level(levelf)
    {
        int size = pow(2, level) + 2*ghost;
        col = size; row = size;
        x.resize((col * row), 0.0);
        y.resize((col * row), 0.0);
        xv.resize((col * row), 0.0);
        yv.resize((col * row), 0.0);

        delta = fmax(width, height) / ((double) pow(2, levelf));
    }

    void resize (double widthf, double heightf, int levelf)
    {
        width = widthf, height = heightf, level = levelf;
        int size = pow(2, level) + 2*ghost;
        col = size; row = size;
        x.resize(row, col, 0.0);
        y.resize(row, col, 0.0);
        xv.resize(row, col, 0.0);
        yv.resize(row, col, 0.0);

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
        allVectorFields.push_back(this);
    }

    void resize(int rows, int cols, int newVal = 0.0)
    {
        nx = cols; ny = rows;
        x.resize(rows, cols, newVal); y.resize(rows, cols, newVal);
    }

    bool empty();

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

    bool empty();

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

    bool empty();

    static std::vector<ScalarField*> allScalarFields;

    std::array<std::function<double()>, NUM_BOUNDARIES> boundary;

    std::array<bool, NUM_BOUNDARIES> dirichlet;
};


struct VertexField
{
    int nx = 0, ny = 0;
    std::vector<double> data;

    VertexField() :
        nx(DCOL), ny(DROW), data(DCOL * DROW, 0.0)
    {
        allVertexFields.push_back(this);
        for (int k = LEFT; k < NUM_BOUNDARIES; ++k)
        {
            dirichlet[k] = false;
        }
    }

    VertexField(int rows, int cols, double initVal = 0.0)
        : nx(grid.col), ny(grid.row), data(grid.col * grid.row, initVal)
    {
        allVertexFields.push_back(this);
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

    bool empty();

    static std::vector<VertexField*> allVertexFields;

    std::array<std::function<double()>, NUM_BOUNDARIES> boundary;

    std::array<bool, NUM_BOUNDARIES> dirichlet;
};



void init_grid (double width, double height, int level);


bool is_ghost_cell (int i, int j);


bool is_real (Point pc);


Point locate (double x, double y);

template <typename T>
void dirichlet_boundary (T& vf, Boundaries boundary, double val)
{
    vf.boundary.dirichlet_x[boundary] = true;
    vf.boundary.dirichlet_y[boundary] = true;
    vf.boundary.x[boundary] = [val]() -> double {return val;};
    vf.boundary.y[boundary] = [val]() -> double {return val;};
}

template <typename T>
void dirichlet_boundary_x (T& vf, Boundaries boundary, double val)
{
    vf.boundary.dirichlet_x[boundary] = true;
    vf.boundary.x[boundary] = [val]() -> double {return val;};
}

template <typename T>
void dirichlet_boundary_y (T& vf, Boundaries boundary, double val)
{
    vf.boundary.dirichlet_y[boundary] = true;
    vf.boundary.y[boundary] = [val]() -> double {return val;};
}


void outlet_boundary (ScalarField& sf, Boundaries boundary);


void update_boundary_impl (VectorField& vf);


void update_boundary_impl (ScalarField& vf);


void update_boundary_impl (FaceVectorField& vf);


template <typename... Fields>
void update_boundary (Fields&... fields) {
    (update_boundary_impl(fields), ...);
}


void fill_facefield (FaceVectorField& uf, VectorField& u);


template <typename T>
double interpolate_bi (Point point, T field, double xp, double yp)
{
    int i = point.i, j = point.j;
    double xi = (xp - grid.x(i,j))/delta;
    double yi = (xp - grid.x(i,j))/delta;

    int offx = sign(xi), offy = sign(yi);
    xi = fabs(xi), yi = fabs(yi);
    return ((field(i,j)*(1. - xi) + field(i + offx,j)*xi)*(1. - yi) +
            (field(i,j+offy)*(1. - xi) + field(i+offx,j+offy)*xi)*yi);
    //return ((field));
}

extern Array2D& x, y; // doesn't work

#include "events.h"
#endif // 2D_GRID_H
