#include "user_funs.h"
#include <cmath>

#define PI 3.14

matrix fun1(matrix x, matrix ud1, matrix ud2)
{
    static uint32_t call_count = 0;
    call_count++;
    matrix y;
    y = -cos(0.1*x()) * exp(-pow((0.1 * x() - 2 * PI),2)) + 0.002*pow((0.1*x()),2);
    return y;
}

matrix fun2(matrix x, matrix ud1, matrix ud2){
    static uint32_t call_count = 0;
    call_count++;
    matrix y;
    matrix y_ref(2, new double[2]{3.14, 0});
    matrix y0(2,1);
    matrix* y = solve_ode(df, 0, 0.1, 100, y0, y_ref, x);
    y=0;
    int n=get_len(y[0]);

    for(int i=0; i<n;i++)
    {
        y=y+10*pow(y_ref(0)-y[1](i,0),2) + pow(y_ref(1)-y[1](i,1),2) + pow(x(0)*(y_ref(0)-y[1](i,0))) + x(1)*(y_ref(1)-y[1](i,1),2);        
    }
    y=y*0.1;
    return y;
}

matrix df(double t, matrix Y, matrix ud1, matrix ud2)
{
    double l = 0.5, mr = 1.0, mc = 9.0,b=0.5;
    double I = (1.0/3.0)*(mr*l*l) + mc*l*l;
    matrix dy(3,1);
    dy(0) = Y(1);
    dy(1) = (ud2(0) + (ud1(0)-Y(0))+ud2(1)+(ud1(1)-Y(1)-b+Y(1)))/I;
    return dy;
}