#pragma once

#include"ode_solver.h"
#include "matrix.h"


matrix fun1(matrix,matrix=NAN,matrix=NAN);
matrix funD(double, matrix,matrix,matrix);
matrix funRP(matrix,matrix=NAN,matrix=NAN);

matrix fun2(matrix,matrix=NAN,matrix=NAN);
matrix df(double,matrix,matrix,matrix);
matrix fun2RP(matrix,matrix=NAN,matrix=NAN);

matrix fun3(matrix,matrix=NAN,matrix=NAN);
matrix df3(double, matrix,matrix, matrix);
matrix fun3RP(matrix,matrix=NAN,matrix=NAN);

matrix fun4(matrix,matrix=NAN,matrix=NAN);
matrix fun4RP(matrix,matrix=NAN,matrix=NAN);
matrix gf(matrix,matrix=NAN,matrix=NAN);
matrix gfRP(matrix,matrix=NAN,matrix=NAN);
matrix hf(matrix,matrix=NAN,matrix=NAN);

matrix fun5(matrix,matrix=NAN,matrix=NAN);
matrix fun5RP(matrix,matrix=NAN,matrix=NAN);