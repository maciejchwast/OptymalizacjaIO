/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
*********************************************/

#include"opt_alg.h"
#include "user_funs.h"
#include <fstream>
#include <iostream>

void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
        std::cout<<lag(fun1,1e-4,1e-2,1e-07,1e-200,1000)<<std::endl;
        /*double* ans = expansion(fun1,69,2,0.98,100);
        std::cout<<ans[0]<<", "<<ans[1]<<std::endl;
        std::cout<<"Fibbonacci:"<<endl;
        std::cout<<fib(fun1,2,20,1).x<<std::endl;
         */
        fstream s;
        s.open("wyniki.txt", ios::out);
        if(!s) throw runtime_error("chuja sie stworzylo a nie plik");
        /*for(int i = 0; i<100; i++)
        {
            double x0 = rand() % 101;
            double d = -8;
            double alpha = 6.3;
            double epsilon = 1e-5;
            double gamma = 1e-200;
            int Nmax = 1000;
            double *p = expansion(fun1, x0, d, alpha, Nmax);
            cout << "X0 = " << x0<<endl;
            cout << p[0] << "\t" << p[1] << endl;

            solution ep;
            s<<"nr"<<i<<"\t\t"<<x0<<" "<<p[0]<<"\t"<<p[1]<<"\t"<<solution::f_calls;
            s<<"\t";

            string c = "local";
            solution::clear_calls();
            solution opt_F = fib(fun1,p[0],p[1], epsilon);
            s<<"Fibbonacci"<<"\n";
            if(opt_F.x>50){
                c = "global";
            }
            s<<opt_F.x<<" "<<opt_F.y<<"\t"<<solution::f_calls<<"\t"<<c<<"\t";
            */
            /*
            solution::clear_calls();
            c = "local";
            solution opt_L = lag(fun1, p[0], p[1],epsilon,gamma,Nmax);
            s<<"Lagrange"<<"\n";
            if(opt_F.x>50){
                c = "global";
            }
            s<<opt_L.x<<" "<<opt_L.y<<"\t"<<opt_L.f_calls<<"\t"<<c<<endl;
             */
        //}
        //s.close();
    }
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab1()
{

}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
