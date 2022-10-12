#pragma once

#include"ode_solver.h"
#include "matrix.h"


matrix fun1(matrix, matrix=0, matrix=0);

double *
expansion(matrix (*ff)(matrix, matrix, matrix), double x, double d, double alpha, int N_max, matrix ud1, matrix ud2);
