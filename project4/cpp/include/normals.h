#ifndef NORMALS_H
#define NORMALS_H

#include "2D-grid.h"

enum HeightType
{
    COLUMN,
    ROW
};

Vector youngs_normal (Point pc, ScalarField& c);

double col_height (Point pc, ScalarField&c, double height = 3);

double row_height (Point pc, ScalarField&c, double height = 3);

double height_function (Point pc, ScalarField&c, double height = 3, HeightType type = COLUMN);

Vector height_normal (Point pc, ScalarField& c, double height = 3);

Vector mycs_normal (Point pc, ScalarField& c, double height = 3);

#endif