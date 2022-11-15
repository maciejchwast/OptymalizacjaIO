#include "user_funs.h"
#include <cmath>

#define PI 3.14

matrix fun1(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    y = -cos(0.1*x()) * exp(-pow((0.1 * x() - 2 * PI),2)) + 0.002*pow((0.1*x()),2);
    return y;
}
/*
double*
expansion(matrix (*ff)(matrix, matrix, matrix), double x, double d, double alpha, int N_max, matrix ud1, matrix ud2)
{
    static uint32_t call_count = 0;
    call_count++;
    int i = 0;
    double* retval = new double [2];
    double x1 = x + d;
    if(ff(x1,0,0) == ff(x,0,0)){
        retval[0] = x;
        retval[1] = x1;
        return retval;
    }

    if(ff(x1,0,0) >ff(x,0,0)){
        d = -d;
        x1 = x +d;
        if(ff(x1,0,0) >= ff(x,0,0)){
            retval[0] = x1;
            retval[1] = x -d;
            return retval;
        }
    }
    matrix next_x = ff(x,0,0);
    matrix prev_x = next_x;

    do
    {
        if(call_count>N_max){
            throw;
        }
        i = i+1;
        prev_x = next_x;
        next_x = x + pow(alpha,i)*d;

    }
    while(ff(prev_x,0,0) <= ff(next_x, 0, 0));

    if(d>0)
    {
        retval[0] = det(prev_x);
        retval[1] = det(next_x);
        return retval;
    }
    retval[0] = det(next_x);
    retval[1] = det(prev_x);
    return retval;
    }*/
