#pragma once

#include"ode_solver.h"
#include "matrix.h"


matrix fun1(matrix, matrix, matrix);
matrix funD(double, matrix,matrix,matrix);
matrix funRP(matrix, matrix, matrix);

matrix fun2(matrix,matrix,matrix);
matrix df(double,matrix,matrix,matrix);
matrix fun2RP(matrix, matrix, matrix);

