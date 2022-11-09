#include"opt_alg.h"

double *
expansion(matrix(*ff)(matrix, matrix, matrix), double x, double d, double alpha, int N_max, matrix ud1, matrix ud2) {
    try {
        //double* p = new double[2]{ 0,0 };
        //Tu wpisz kod funkcji
        static uint32_t call_count = 0;
        call_count++;
        int i = 0;
        double *p = new double[2];
        double x1 = x + d;
        if (ff(x1, 0, 0) == ff(x, 0, 0)) {
            p[0] = x;
            p[1] = x1;
            return p;
        }

        if (ff(x1, 0, 0) > ff(x, 0, 0)) {
            d = -d;
            x1 = x + d;
            if (ff(x1, 0, 0) >= ff(x, 0, 0)) {
                p[0] = x1;
                p[1] = x - d;
                return p;
            }
        }
        matrix next_x = ff(x, 0, 0);
        matrix prev_x = next_x;

        do {
            if (call_count > N_max) {
                throw;
            }
            i = i + 1;
            prev_x = next_x;
            next_x = x + pow(alpha, i) * d;

        } while (ff(prev_x, 0, 0) <= ff(next_x, 0, 0));

        if (d > 0) {
            p[0] = det(prev_x);
            p[1] = det(next_x);
            return p;
        }
        p[0] = det(next_x);
        p[1] = det(prev_x);
        std::cout << p[0] << std::endl << p[1] << std::endl;
        return p;
        //return p;
    }
    catch (string ex_info) {
        throw ("double* expansion(...):\n" + ex_info);
    }
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji
        double next_c;
        int fi[100];
        int k=100;
        fi[0]=1;
        fi[1]=1;
        for(int i=2;i<100;i++)
        {
            fi[i]=fi[i-1]+fi[i-2];
        }

        double c = b - fi[k-1]/(fi[k]*b-a);
        double d = a + b - c;
        double next_a;
        double next_b;
        for(int i =0; i<k-3;i++)
        {
            if(ff(c,0,0)<ff(d,0,0))
            {
                next_a = a;
                next_b = d;
            }
            else {
                next_b = b;
                next_a = a;
            }

           next_c = next_b-fi[k-i-2]/fi[k-i-1]*(next_b-next_a);
        }
        std::cout<<"fibo:"<<std::endl;
        for(int i =0; i<k-3;i++)
        {
            if(ff(c,0,0)<ff(d,0,0))
            {
                next_a = a;
                next_b = d;
            }
            else {
                next_b = b;
                next_a = a;
            }

            next_c = next_b-fi[k-i-2]/fi[k-i-1]*(next_b-next_a);
            if (fi[i] > (b - a) / epsilon) {
                std::cout << i << std::endl;
                break;
            }
        }
        double z;
        //return z*=next_c;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution
lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1,
    matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution lag(...):\n" + ex_info);
    }
}

solution
HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1,
   matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji
        solution X, XB, XB_old;
        bool failed = false;
        do {
            X.fit_fun(ff, ud1, ud2);
            XB.fit_fun(ff, ud1, ud2);
            if (X.y < XB.y) {
                do {
                    XB_old.x = XB.x;
                    XB.x = X.x;
                    X.x = 2 * XB.x - XB_old.x;
                    X = HJ_trial(ff, XB, s);
                    if (X.f_calls > Nmax) {
                        failed = true;
                        throw std::runtime_error("Calls exceeded!");
                    }
                } while (X.y >= XB.y);
            } else {
                s = alpha * s;
            }
            if (X.f_calls > Nmax) {
                failed = true;
                throw std::runtime_error("Calls exceeded!");
            }
        } while (s < epsilon);
        X.fit_fun(ff, ud1, ud2);
        Xopt = X;
        if (failed) Xopt.flag = 1;
        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution HJ(...):\n" + ex_info);
    }
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2) {
    try {
        //Tu wpisz kod funkcji
        int n = get_dim(XB);
        matrix E = ident_mat(n);
        solution X;
        for (int j = 0; j < n; j++) {
            X.x = XB.x + s * E[j];
            X.fit_fun(ff, ud1, ud2);
            if (X.y < XB.y) {
                XB = X;
            } else {
                X.x = XB.x = s * E[j];
                X.fit_fun(ff, ud1, ud2);
                if (X.y < XB.y) {
                    XB = X;
                }
            }
        }
        return XB;
    }
    catch (string ex_info) {
        throw ("solution HJ_trial(...):\n" + ex_info);
    }
}

solution
Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax,
      matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji
        solution XB(x0), X;
        int n = get_dim(XB);
        matrix lambda(n, 1), p(n, 1), s(s0), D = ident_mat(n);


        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution Rosen(...):\n" + ex_info);
    }
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1,
             matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution pen(...):\n" + ex_info);
    }
}

solution
sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta,
       double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution sym_NM(...):\n" + ex_info);
    }
}

solution
SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon,
   int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution SD(...):\n" + ex_info);
    }
}

solution
CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon,
   int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution CG(...):\n" + ex_info);
    }
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
                matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1,
                matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution Newton(...):\n" + ex_info);
    }
}

solution
golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution golden(...):\n" + ex_info);
    }
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution Powell(...):\n" + ex_info);
    }
}

solution
EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix limits, int mi, int lambda, matrix sigma0, double epsilon,
   int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution EA(...):\n" + ex_info);
    }
}
