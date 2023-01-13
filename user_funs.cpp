#include "user_funs.h"
#include <cmath>
#include "solution.h"

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
    //Y[0] to czas Y[1] to rozwiazanie
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

matrix fun3(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    double a = PI * sqrt((pow(x(0)/PI, 2)+pow(x(1)/PI,2)));
    y = sin(a)/a;
    if(ud2(1)>1) //zewnetrzna funkcja kary
    {
        if(-x(0) + 1>0)
        {
            y = y + ud2(0) * pow(-x(0) + 1, 2);   
        }
        if(-x(1) + 1>0)
        {
            y = y + ud2(0) * pow(-x(1) + 1, 2);
        }
        if(norm(x) - ud1(0) > 0)
        {
            y = y + ud2(0) * pow(norm(x) - ud1(0), 2);
        }
    }
    else //wewnetrzna funkcja kary
    {
        if(-x(0) + 1 > 0)
        {
            y = 1e10;
        }else {
            y = y - ud2(0) / (-x(0) + 1);
        }
        if(-x(1) + 1 > 0)
        {
            y = 1e10;   
        } else {
            y = y + ud2(0) / (-x(1) + 1);
        }
        if(norm(x) - ud1(0) > 0)
        {
            y = 1e10;   
        } else {
            y = y + ud2(0) / (norm(x) - ud1(0));   
        } 
    }
    return y;
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2)
{
    double c = 0.47, r = 0.12, m = 0.6, ro = 1.2, g = 9.81;
    double s = 3.14 * r * r;
    double omega = (ud1)(0);
    double Dx = 0.5 * c * ro * s * abs(Y(1))*Y(1);
    double Dy = 0.5 * c * ro * s * abs(Y(3))*Y(3);
    double FMx = 3.14 * ro * Y(3) * omega * pow(r,3);
    double FMy = 3.14 * ro * Y(1) * omega * pow(r,3);

    matrix dY(4,1);
    dY(0) = Y(1);
    dY(1) = (-Dx -FMx)/m;
    dY(2) = Y(3);
    dY(3) = (-m*g -Dy - FMy)/m;
    return dY;
}
                         
matrix fun3RP(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    matrix Y0(4,new double[4]{0,x(0),100,0});
    matrix omega = x(1);
    matrix *Y = solve_ode(df3,0,0.01,7,Y0,omega);
    int n = get_len(Y[0]);
    int i0 = 0, i50 =0;

    for (int i = 0; i < n; ++i) {
        if(abs(Y[1](i,2) - 50)<abs(Y[1](i50,2)-50)){
            i50 = i;
        }
        if(abs(Y[1](i,2))<abs(Y[1](i0,2))){
            i0=i;
        }
    }

    y = - Y[1](i0,0);
    // 3 ograniczenia
    if(abs(x(0)) - 10>0){
        y = y + (ud2)(0)*pow(abs(x(0))-10,2);
    }
    if (abs(x(1)) - 20>0){
        y = y+(ud2)(0)* pow(abs(x(1))-20,2);
    }
    if(abs(Y[1](i50,0)-5)-1>0){
        y = y+(ud2)(0)* pow(abs(Y[1](i50,0)-5)-1,2);
    }
    return y;
}

matrix fun4(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    if(isnan(ud2(0,0)))
    {
        y=pow(x(0)+2*x(1)-7,2)+ pow(2*x(0)+x(1) - 5,2);
    }
    else
    {
        y=fun4(ud2(0)+x*ud2(1),ud1,ud2);
    }
    return y;
}

matrix fun4RP(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    int m = 100, n = get_len(x);
    static matrix X(n,m), Y(1,m);
    if(solution::f_calls == 1)
    {
        ifstream S("XData.txt");
        S >> X;
        S.close();
        S.open("YData.txt");
        S >> Y;
        S.close();
    }
    double h;
    y = 0;

    for(int i = 0; i < m ; i++){
        h=(trans(x)*X[i])();
        h = 1 /(1+exp(-h)); //hipoteza
        y = y - Y(0,i) * log(h) - (1-Y(0,i)) * log(1-h);
    }
    y = y/m;
}

matrix gf(matrix x,matrix ud1,matrix ud2)
{
    matrix g(2,1);
    g(0) = 10*x(0) + 8*x(1) - 34;
    g(1) = 8*x(0) +10*x(1)-38;
    return g;
}

matrix hf(matrix x, matrix ud1, matrix ud2)
{
    matrix H(2,2); // hesjan jest sta³y funkcja nie zmienia swojej wypukloœci ma tylko jedo minimum
    H(0,0) = H(1,1) = 10;
    H(0,1)=H(1,0) = 8;
    return H;
}

matrix gfRP(matrix x, matrix ud1, matrix ud2)
{
    int m = 100, n = get_len(x);
    static matrix X(n,m), Y(1,m);
    if(solution::g_calls == 1)
    {   ifstream S("XData.txt");
        S >> X;
        S.close();
        cout << X << endl;
        S.open("YData.txt");
        S >> Y;
        S.close();
        cout << Y << endl;

    }

    double h;
    matrix g(n,1);

    for(int j = 0; j < n; j++)
    {
        for(int i = 0; i < m; i++)
        {
            h = (trans(x)*X[i])();
            h = 1/(1+ exp(-h));
            g(j) = g(j) + X(j,i) * (h- Y(0,i));
        }
        g(j) = g(j)/m;
    }
}

matrix fun5(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    if(ud2 == NULL)
    {
        y = matrix(2,1); //2 funkcje celu
        y(0) = ud1[0]() * pow(x(0)-2, 2) + pow(x(1) -2,2);
        y(1) = 1/ud1[0]() * (pow(x(0) + 2,2) + pow(x(1) +2,2));
    }
    else // y = g(alfa)
    {
        solution T;
        T.x = ud2[0] +x*ud2[1]; //ud2[0] - punkt w ktorym sie obecnie znajdujemy
        T.fit_fun(fun5,ud1,ud2);
        y = ud1[1]()*T.y(0) + (1 - ud1[1]()) * T.y(1);
        solution::f_calls--;
    }

}

matrix fun5RP(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    if(ud2 == NULL)
    {
        y = matrix(3,1);
        double ro = 7800, p = 1e3, E = 207e9;
        y(0) = ro *x(0) *3.14 * pow(x(1),2)/4; //masa
        y(1) = 64 * p * pow(x(0),3)/(3*E*3.14*pow(x(1),4)); //ugiecie
        y(2) = 32 * p *x(0)/(3.14*pow(x(1),3)); //naprezenie
    }
    else
    {
        matrix yt, xt = ud2[0] + x * ud2[1];
        yt = fun5RP(xt,ud1,NULL);
        y = ud1 * (yt(0) - 0.06)/(1.53 - 0.06) + (1 -ud1) * (yt(1) - 5.25e-6)/(0.0032 - 5.25e-6);

        double c =1e10;
        if(xt(0)<0.2)
        {
            y=y+c*pow(0.2 -xt(0),2);
        } 
        if(xt(0)>1)
        {
            y=y+c*pow(xt(0)-1,2);
        } 
        if(xt(0)<0.01)
        {
            y=y+c*pow(0.01-xt(1),2);
        } 
        if(xt(0)>0.05)
        {
            y=y+c*pow(xt(1)-0.05,2);
        } 
        if(xt(0)>0.005)
        {
            y=y+c*pow(yt(1)-0.005,2);
        } 
        if(xt(0)<300e6)
        {
            y=y+c*pow(yt(2)-300e6,2);
        } 
    }
}