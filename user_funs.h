#pragma once

#include"ode_solver.h"
#include "matrix.h"


matrix fun1(matrix, matrix, matrix);
matrix funD(double, matrix,matrix,matrix);
matrix funRP(matrix, matrix, matrix);

matrix fun2(matrix,matrix,matrix);
matrix df(double,matrix,matrix,matrix);
matrix fun2RP(matrix, matrix, matrix);

matrix fun3(matrix,matrix,matrix);
matrix df3(double, matrix,matrix, matrix);
matrix fun3RP(matrix,matrix,matrix);


matrix fun5(matrix, matrix, matrix);
matrix fun5RP(matrix, matrix, matrix);