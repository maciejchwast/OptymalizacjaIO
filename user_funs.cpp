#include "user_funs.h"
#include <cmath>

#define PI 3.14

matrix fun1(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    y = -cos(0.1*x()) * exp(-pow((0.1 * x() - 2 * PI),2)) + 0.002*pow((0.1*x()),2);
    return y;
}


matrix funRP(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    matrix X0 = matrix(3, new double [3] {5,1,10});
    matrix* X = solve_ode(funD, 0, 1, 1000, X0, ud1,x);
    int l = get_len(X[0]);
    double max = X[1](0,2);
    for(int i =0; i<l; i++)
    {
        if(max < X[1](i,2)) max = X[1](i,2);
    }
    y = abs(max-50);
    return y;
}
matrix funD(double t, matrix X, matrix ud1, matrix ud2)
{
    double Pa=0.75, Va = 5, Ta = 90; //dane zbiornika A
    double Pb = 1, Vb = 1, Tb = 10, Db = 0.003656652; //dane zbiornika B
    double Tin = 10, Fin = 0.01; //dane wlewajacej sie wody F- predkosc przelewu 0.01 bo zmiana z dm3 na m3
    //dane do wzoru na zmiane objetosci
    double a =0.98, b = 0.63, g = 9.81;

    matrix Dx(3,1);
    double Fa, Fb;
    if(X(0)>0)
    {
        Fa = (-1.0* a) * b * m2d(ud2) * sqrt(2 * g * X(0)/Pa);
    }
    else Fa = 0;
    if(X(1)>0)
    {
        Fb = (-1.0*a)* b * Db * sqrt(2 * g * X(1)/Pb);
    }
    else Fb = 0;
    double Tzm = Fin/X(1) * (Tin - X(2)) -1.0 * Fa / X(1) * (Ta - X(2)); //zmiana temp wody

    Dx(0) = Fa;
    Dx(1) = Fb - Fa + Fin;
    Dx(2) = Tzm;

    return Dx;
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
