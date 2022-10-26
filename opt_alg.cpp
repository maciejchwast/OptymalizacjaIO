#include"opt_alg.h"
#include "user_funs.h"
double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alfa, int Nmax, matrix ud1, matrix ud2) {
    int i = 0;
    double x1 = x0 + d;

    if (ff(x1, NULL, NULL) == ff(x0, NULL, NULL)) {
        double* tab = new double[2];
        tab[0] = x0;
        tab[1] = x1;
        return tab;
    }

    if (ff(x1, NULL, NULL) > ff(x0, NULL, NULL)) {
        d = -d;
        x1 = x0 + d;
        if (ff(x1, NULL, NULL) >= ff(x0, NULL, NULL)) {
            double* tab = new double[2];
            tab[0] = x1;
            tab[1] = x0 - d;
            //cout << "petla: " << ff(x0, NULL, NULL);
            //cout << "petla: " << ff(x1, NULL, NULL);
            return tab;
        }
    }

    vector<double> vec;
    vec.push_back(x0);
    vec.push_back(x1);

    do {
        if (i > Nmax) {
            double* tab = new double;
            *tab = 404;
            return tab;
        }

        i++;

        vec.push_back(vec.at(vec.size() - 1) + alfa * d);

    } while (ff(vec.at(vec.size() - 2), NULL, NULL) >= ff(vec.at(vec.size() - 1), NULL, NULL));

    if (d > 0) {
        double* tab = new double[2];
        tab[0] = vec.at(vec.size() - 3);
        tab[1] = vec.at(vec.size() - 1);
        return tab;

    }
    else {
        double* tab = new double[2];
        tab[1] = vec.at(vec.size() - 1);
        tab[0] = vec.at(vec.size() - 3);
        return tab;
    }
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon) {
    vector<double> fibNumbers;
    fibNumbers.push_back(0);
    fibNumbers.push_back(1);
    for (int i = 2; i < 100; i++) {
        fibNumbers.push_back(fibNumbers.at(i - 1) + fibNumbers.at(i - 2));
    }

    int k = 3;
    for (int i = 0; i < fibNumbers.size(); i++) {
        if (fibNumbers.at(i) > (b - a) / epsilon) {
            k = i;
            break;
        }
    }

    // cout << k << endl; // jest ok

    double a1, b1, c, d;
    a1 = a;
    b1 = b;
    c = b1 - fibNumbers.at(k - 1) / fibNumbers.at(k) * (b1 - a1);
    d = a + b1 - c;

    for (int i = 0; i <= k - 3; i++) {
        if (fun(c, 0, 0) < fun(d, 0, 0)) {
            b1 = d;
        }
        else {
            a1 = c;
        }
        c = b1 - fibNumbers.at(k - i - 2) / fibNumbers.at(k - i - 1) * (b1 - a1);
        d = a1 + b1 - c;
    }

    return c;
}


solution lag(matrix(*ff)(matrix, matrix, matrix), double aInput, double bInput, double cInput, double eps, double gamma, int Nmax) {
    int i = 0;
    double a = aInput;
    double b = bInput;
    double c = cInput;
    double l, m;
    double fcalls = 0;
    double d0, d=2;
    do {

        d0 = d;
        fcalls++;
        l = m2d(fun(a,0,0) * (b * b - c * c) + fun(b,0,0) * (c * c - a * a) + fun(c,0,0) * (a * a - b * b));
        m = m2d(fun(a,0,0) * (b - c) + fun(b,0,0) * (c - a) + fun(c,0,0)*(a - b));

        if (m <= 0) return 404;

        d = 0.5 * l / m;

        if (a<d && c>d) {
            if (fun(d) < fun(c)) {
                a = a;
                b = c;
                c = d;
            }
            else {
                a = d;
                c = c;
                b = b;
            }
        }
        else {
            if (c<d && b>d)
                if (fun(d) < fun(c)) {
                    a = c;
                    c = d;
                    b = b;
                }
                else {
                    a = a;
                    c = c;
                    b = d;
                }
            else return 4042;
        }
        i++;
        if (fcalls > Nmax)
            return 40422;

    } //while (b - a < eps || abs(d - d0) < gamma);
    while (b - a >= eps && abs(d - d0) >= gamma);
    return d;
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix limits, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
