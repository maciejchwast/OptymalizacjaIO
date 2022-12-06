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

matrix fun2(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    y=pow(x(0),2)+pow(x(1),2)- cos(2.5*3.14*x(0)) - cos(2.5 * 3.14 *x(1)) +2;
    return y;
}

matrix df(double t, matrix Y, matrix ud1, matrix ud2)
{
    double mr =1, mc =9, l = 0.5, b=0.5, a_ref = 3.14, o_ref = 0;
    double I = (mr * l * l)/3 + mc * l * l, k1 = (ud2)(0), k2 = (ud2)(1);
    double M = k1 * (a_ref -Y(0)) + k2*(o_ref - Y(1));

    matrix dY(2,1);
    dY(0) = Y(1);
    dY(1) = (M - b*Y(1))/I;
    return dY;
}

matrix fun2RP(matrix x, matrix ud1, matrix ud2)
{
    //Y[0] to jest czas Y[1] to jest rozwiązanie
    matrix y;
    matrix Y0(2,1);
    matrix *Y = solve_ode(df,0,0.1,100,Y0,ud1,x);
    double a_ref = 3.14 , o_ref = 0;
    int n = get_len(Y[0]);
    y = 0;
    for (int i = 0; i < n; ++i) {
        y = y+10* pow(a_ref-Y[1](i,0), 2)+pow(o_ref - Y[1](i,1),2)+ pow(x(0)*(a_ref-Y[1](i,0))+ x(1)*(o_ref - Y[1](i,1)),2);
    }
    y = y*0.1;
    return y;
}
