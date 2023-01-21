/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
*********************************************/

#include"opt_alg.h"
#include "user_funs.h"
#include <fstream>
#include <iostream>

void proj1();
void proj2();
void proj3();
void proj3RP();
void proj4();
void proj4RP();
void proj5();
void proj5RP();


int main()
{
	try
    {
        proj5();
    }

	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}

	return 0;
}

void proj1()
{
    std::cout<<"Expansion"<<std::endl;
    double* ans = expansion(fun1,69,2,0.98,100);
    std::cout<<ans[0]<<", "<<ans[1]<<std::endl<<std::endl;
    std::cout<<"Fibbonacci:"<<endl;
    std::cout<<fib(fun1,-100,100,1)<<std::endl;
    solution::clear_calls();
    std::cout<<"Lagrange:"<<std::endl;
    std::cout<<lag(fun1,-12,34,0.001,1e-6,1000)<<std::endl;
    solution::clear_calls();
    std::cout<<"Rzeczywisty:"<<std::endl;

    std::cout<<fib(funRP,1e-4,1e-2,1e-7)<<std::endl;
    solution::clear_calls();
    std::cout<<lag(funRP,1e-4,1e-2,1e-7,1e-200,1000)<<std::endl;

    fstream s;
    s.open("wyniki.txt", ios::out);
    if(!s) throw runtime_error("nic sie nie stworzylo");
    for(int i = 0; i<100; i++)
    {
        double x0 = 13.21 +(3.0/7.0)*i;
        double d = -8;
        double alpha = 6.3;
        double epsilon = 1e-5;
        double gamma = 1e-200;
        int Nmax = 1000;

        double *p = expansion(fun1, x0, d, alpha, Nmax);
        //cout << "X0 = " << x0<<endl;
        //cout << p[0] << "\t" << p[1] << endl;
        solution ep;
        //s<<"nr"<<i<<"\t\t"<<x0<<" "<<p[0]<<"\t"<<p[1]<<"\t"<<ep.f_calls;
        //s<<"\t";
        string c = "local";

        solution::clear_calls();
        solution opt_F = fib(fun1,-100,100, epsilon);
        //s<<"Fibbonacci"<<"\n";
        if(opt_F.x>50){
            c = "global";
        }
        //s<<opt_F.x<<" "<<opt_F.y<<"\t"<<solution::f_calls<<"\t"<<c<<"\t";
        //s<<opt_F.ud;

        solution::clear_calls();
        c = "local";
        solution opt_L = lag(fun1, -100, 100,epsilon,gamma,Nmax);
        s<<"Lagrange"<<"\n";
        if(opt_L.x>50){
            c = "global";
        }
        s<<opt_L.x<<" "<<opt_L.y<<"\t"<<opt_L.f_calls<<"\t"<<c<<endl;
        //s<<opt_L.ud;
    }
    s.close();


    fstream g;
    g.open("ob.txt", ios::out);
    double x0 = 0.0001;
    double d = 0.001;
    double alpha = 1.5;
    double epsilon = 1e-7;
    double gamma = 1e-200;
    int Nmax = 1000;

    double *p = expansion(funRP, x0, d, alpha, Nmax);

    solution::clear_calls();
    solution opt_F = fib(funRP,p[0], p[1], epsilon);
    //cout << opt_F << endl;
    solution::clear_calls();
    solution opt_L = lag(funRP, p[0], p[1], epsilon, gamma, Nmax);
    //g << opt_L;
    cout << opt_L << endl;
    //cout << endl;
    matrix X0 = matrix(3, new double[3]{ 5,1,10 });
    matrix* X = solve_ode(funD,0, 1, 1000, X0,NULL,opt_L.x);
    g<<X[1];
}

void proj2() {
    /*double s = 0.0012312312, alpha_HJ = 0.5, alpha_R = 2, beta = 0.5, epsilon = 1e-3;
    //s - dlugosc kroku
    int nMax = 1000;
    matrix x0, s0(2, 1, s);
    solution opt_HJ;
    solution opt_R;
    string k = "NIE";
    for (int i = 0; i < 100; ++i) {
        solution::clear_calls();
        x0 = 2 * rand_mat(2, 1) - 1;

        opt_HJ = HJ(fun2, x0, s, alpha_HJ, epsilon, nMax);

        if (opt_HJ.y<0.01 && opt_HJ.y>(-0.01))k = "TAK"; // warunek miejsca zerowego
        cout << x0(0) << " " << x0(1) << " " << opt_HJ.x(0) << " " << opt_HJ.x(1) << " " << opt_HJ.y << " "
             << solution::f_calls << " " << k << "\n";
        k = "NIE";

        solution::clear_calls();
        opt_R = Rosen(fun2, x0, s0, alpha_R, beta, epsilon, nMax);
        if (opt_R.y<0.01 && opt_R.y>(-0.01))k = "TAK"; // warunek miejsca zerowego
        cout << opt_R.x(0) << " " << opt_R.x(1) << " " << opt_R.y << " " << solution::f_calls << " " << k << endl;
        //cout << opt_R<<endl;
        k = "NIE";
    }
    double s = 0.1379, alpha_HJ = 0.5, alpha_R = 2, beta = 0.5, epsilon = 1e-3;
    int nMax = 1000;

    matrix x0 = 2 * rand_mat(2, 1) - 1, s0(2, 1, s);
    matrix XS_HJ = trans(x0);
    solution::clear_calls();
    solution opt_HJ = HJ(fun2,x0, s, alpha_HJ, epsilon, nMax, XS_HJ);

    cout<<"H:"<<endl;
    cout<<XS_HJ[0]<<"\n"<<XS_HJ[1];
    cout<<opt_HJ.f_calls;
    matrix XS_R = trans(x0);
    solution::clear_calls();
    solution opt_R = Rosen(fun2,x0, s0, alpha_R, beta, epsilon, nMax,XS_R);
    cout<<"\nR:"<<endl;
    cout<<XS_R[0]<<"\n\n"<<XS_R[1];
    cout<<opt_R.f_calls;*/
    matrix x0(2,1,1);
    double s = 0.1379, alpha_HJ = 0.2, alpha_R = 2, beta = 0.5, epsilon = 1e-3;
    int nMax = 1000;
    matrix s0(2, 1, s);

    solution::clear_calls();
    solution opt_HJ = HJ(fun2RP,x0, s, alpha_HJ, epsilon, nMax);
    //cout<<opt_HJ<<endl;

    solution::clear_calls();
    solution opt_R = Rosen(fun2RP,x0, s0, alpha_R, beta, epsilon, nMax);
    //cout<<opt_R<<endl;

    matrix start = matrix(2, new double[2]{0,0});
    matrix* w_hj = solve_ode(df,0,0.1,100,start,NULL,opt_HJ.x);
    //cout<<w_hj[1]<<endl;

    matrix* w_R = solve_ode(df,0,0.1,100,start,NULL,opt_R.x);
    cout<<w_R[1]<<endl;

}

void proj3()
{
    //dc dla wewnętrznej funkcji kary 0,5
    //dc dla zewnętrznej funkcji kary 2

    double c0 = 9, dc = 2, epsilon = 1e-4;
    int Nmax = 5000;
    matrix x0, a=4; //macierz a

    do{
        x0 = 4* rand_mat(2,1)+1;
    }while(norm(x0)>a);
    cout<<x0<<endl<<endl;

    solution opt_zew = pen(fun3,x0,c0,dc,epsilon,Nmax,a);
    cout<<opt_zew<<endl;
    cout<<norm(opt_zew.x)<<endl; // w pliku xls kolumna r to odleglosc
    // 4,00001 wskazuje na to że jest poza przedziałem

    //-----------------------------------------------------------------

    for (int i = 0; i < 100; ++i) {
        //losowanie punktu startowego
        do{
            x0 = 4* rand_mat(2,1)+1;
        }while(norm(x0)>a);

        //zewnetrzne
        dc = 2;
        solution::clear_calls();
        solution opt_zewn = pen(fun3,x0,c0,dc,epsilon,Nmax,a);
        double r_zewn = norm(opt_zewn.x);
        cout<<x0(0)<<"\t"<<x0(1)<<"\t"<<opt_zewn.x(0)<<"\t"<<opt_zewn.x(1)<<"\t"<<r_zewn<<"\t"<<opt_zewn.y<<"\t"<<solution::f_calls<<"\t";

        //wewnetrzne
        dc = 0.5;
        solution::clear_calls();
        solution opt_wewn = pen(fun3,x0,c0,dc,epsilon,Nmax,a);
        double r_wewn = norm(opt_wewn.x);
        cout<<opt_wewn.x(0)<<"\t"<<opt_wewn.x(1)<<"\t"<<r_wewn<<"\t"<<opt_zewn.y<<"\t"<<solution::f_calls<<"\n";
    }
}

void proj3RP()
{
    matrix x(2,1,1),c=0.47, a =4;
    double c0 = 9, dc = 2, epsilon = 1e-5;
    int Nmax = 10000;
    solution opt_zewn = pen(fun3RP,x,c0,dc,epsilon,Nmax,a);
    opt_zewn.fit_fun(fun3RP,NULL, c);
    cout<<opt_zewn<<endl;

    cout<<opt_zewn.x(1)<<endl;

    matrix Y(4, new double[4]{0,opt_zewn.x(0),100,0});
    matrix ud(opt_zewn.x(1));
    matrix * R = solve_ode(df3,0,0.01,7,Y,ud);
    cout<<R[1];


}

void proj4()
{
    fstream s;
    s.open("lab4.txt", ios::out);
    double h0 = -1, epsilon = 1e-5;
    int Nmax = 10000;
    solution optSD, optCG, optN;
    for (int i = 0; i < 100; ++i) {
        matrix x0 = 20* rand_mat(2,1)-10;
        cout<<x0(0)<<"\t"<<x0(1)<<"\t";
        optSD = SD(fun4,gf,x0,h0,epsilon,Nmax);
        cout<<optSD.x(0)<<"\t"<<optSD.x(1)<<"\t"<<optSD.y(0)<<"\t"<<solution::f_calls<<"\t"<<solution::g_calls<<"\t";
        solution::clear_calls();

        optCG = CG(fun4,gf,x0,h0,epsilon,Nmax);
        cout<<optCG.x(0)<<"\t"<<optCG.x(1)<<"\t"<<optCG.y(0)<<"\t"<<solution::f_calls<<"\t"<<solution::g_calls<<"\t";
        solution::clear_calls();

        optN = Newton(fun4,gf,hf,x0,h0,epsilon,Nmax);
        cout<<optN.x(0)<<"\t"<<optN.x(1)<<"\t"<<optN.y(0)<<"\t"<<solution::f_calls<<"\t"<<solution::g_calls<<"\t";
        cout<<optN.H_calls<<"\n";
        solution::clear_calls();
    }
    /*
    matrix x0 = 20 * rand_mat(2, 1) - 10;
    double h0[] = {0.05, 0.12, -1};
    double epsilon = 1e-3;
    int Nmax = 1000;
    matrix out;
    // Najszybszy spadek
    out = trans(x0);
    solution::clear_calls();
    SD(fun4, gf, x0, h0[0], epsilon, Nmax, out);
    cout << x0<<endl;


    out = trans(x0);
    solution::clear_calls();
    SD(fun4,gf,x0, h0[1], epsilon, Nmax, out);


    out = trans(x0);
    solution::clear_calls();
    SD(fun4,gf,x0, h0[2], epsilon, Nmax, out);

    // Gradient sprzezony
    out = trans(x0);
    solution::clear_calls();
    CG(fun4,gf,x0, h0[0], epsilon, Nmax, out);


    out = trans(x0);
    solution::clear_calls();
    CG(fun4,gf,x0, h0[1], epsilon, Nmax, out);


    out = trans(x0);
    solution::clear_calls();
    CG(fun4,gf,x0, h0[2], epsilon, Nmax, out);


    // Newton
    out = trans(x0);
    solution::clear_calls();
    Newton(fun4,gf,hf,x0, h0[0], epsilon, Nmax, out);


    out = trans(x0);
    solution::clear_calls();
    Newton(fun4,gf,hf,x0, h0[1], epsilon, Nmax, out);


    out = trans(x0);
    solution::clear_calls();
    Newton(fun4,gf,hf,x0, h0[2], epsilon, Nmax, out);
 */

}

void proj4RP()
{
    double h = 0.01, epsilon = 1.e-5;
    int Nmax = 10000;
    matrix x(3, new double[3]{ 0,0,0 });
    solution sol = CG(fun4RP,gfRP,x, h, epsilon, Nmax);
    cout << sol << endl;
    matrix X(3, 100), Y(1, 100);
    ifstream S("Xdata.txt");
    S >> X;
    S.close();
    S.open("Ydata.txt");
    S >> Y;
    S.close();
    double h0, p = 0.0;
    for (int i = 0; i < 100; i++) {
        h0 = (trans(sol.x) * X[i])();
        cout << X[i];
        h0 = 1 / (1 + exp(-h0));
        if ((h0 >= 0.5 && Y[i] == 1) || (h0 < 0.5 && Y[i] == 0)) {
            p++;
        }
    }
    p = p / 100;
    cout << "p= " << p << endl;
}

void proj5()
{
    matrix x0 = 20*rand_mat(2,1) -10, ud(2,1);
    double epsilon = 1e-5, w = 0;
    int Nmax = 5000, a =1;
    ud[0] = a;
    ud[1] = w;

    solution opt = Powell(fun5,x0,epsilon,Nmax,ud);
    //cout<<opt<<endl;

    for (int i = 0; i<101; ++i) {
        if (i == 0) { w = 0; } else { w += 0.01; };
        x0 = 20 * rand_mat(2, 1) - 10;

        a = 1;
        ud[0] = a;
        ud[1] = w;

        solution::f_calls = 0;
        opt = Powell(fun5, x0, epsilon, Nmax, ud);
        cout << x0(0) << "\t" << x0(1) << "\t" << opt.x(0) << "\t" << opt.x(1) << "\t" << opt.y(0) << "\t" << opt.y(1)
             << "\t" << solution::f_calls << "\t";

        a = 10;
        ud[0] = a;
        ud[1] = w;
        solution::f_calls = 0;
        opt = Powell(fun5, x0, epsilon, Nmax, ud);
        cout << opt.x(0) << "\t" << opt.x(1) << "\t" << opt.y(0) << "\t" << opt.y(1) << "\t" << solution::f_calls
             << "\t";

        a = 100;
        ud[0] = a;
        ud[1] = w;
        solution::f_calls = 0;
        opt = Powell(fun5, x0, epsilon, Nmax, ud);
        cout << opt.x(0) << "\t" << opt.x(1) << "\t" << opt.y(0) << "\t" << opt.y(1) << "\t" << solution::f_calls
             << "\n";
    }
}

void proj5RP()
{
    matrix x0(2,1), ud (2,1);
    /*
    //min masa max ugiecie
    x0(0) = 0.2;
    x0(1) = 0.01;
    solution test;
    test.x = x0;
    test.fit_fun(fun5RP);
    cout<<test<<endl;

    //max masa min ugiecie
    x0(0) = 0.2;
    x0(1) = 0.05;
    solution test2;
    test2.x = x0;
    test2.fit_fun(fun5RP);
    cout<<test2<<endl;
*/
    double epsilon = 1e-3;
    int Nmax = 5000;
    solution opt;
    double w=0.01;
    int a;

    for (int i = 0; i < 100; ++i) {
        x0 = 20*rand_mat(2,1) -10;
        ud = (1,1,w);
        opt = Powell(fun5RP,x0,epsilon,Nmax,ud);
        cout<<opt;
        solution::clear_calls();
        w+=0.01;
    }
}
