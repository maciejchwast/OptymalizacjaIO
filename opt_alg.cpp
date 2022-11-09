#include"opt_alg.h"

double * expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int N_max, matrix ud1,
                   matrix ud2) {
    try {
        double* p = new double[2]{ 0,0 };
        //Tu wpisz kod funkcji
        int i = 0;
        double x1  = x0 + d;
        solution X0(x0), X1(x1), X_next(0), X_prev;
        X0.fit_fun(ff,ud1,ud2);
        X1.fit_fun(ff,ud1,ud2);
        if(X0.y == X1.y)
        {
            p[0] = det(X0.x);
            p[1] = det(X1.x);
            return p;
        }
        if(X1.y > X0.y)
        {
            d = -d;
            X1.x = X0.x + d;
            if(X1.y >= X0.y)
            {
                p[0] = det(X1.x);
                p[1] = det(X0.x)-d;
                return p;
            }
        }
        do
        {
            if(X0.f_calls>N_max)
            {
                throw std::runtime_error("Calls number exceeded!");
            }
            X_prev = X_next;
            i++;
            X_next.x = X0.x + pow(alpha, i)*d;
            X_next.fit_fun(ff,ud1,ud2);
        }
        while(X_prev.y <= X_next.y);
        if(d>0)
        {
            p[0] = det(X_prev.x);
            p[1] = det(X_next.x)-d;
            return p;
        }
        p[0] = det(X_next.x);
        p[1] = det(X_prev.x)-d;
        return p;
    }
    catch (string ex_info) {
        throw ("double* expansion(...):\n" + ex_info);
    }
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
        int n = 100;
        int* fib = new int[n] {1, 1};
        for(int i = 2; i < n; i++)
        {
            fib[i] = fib[i-2] + fib[i-1];
        }
        solution A(a), B(b), C, D;
        C.x = B.x - 1.0 * fib[n-2]/fib[n-1] * (B.x - A.x);
        D.x = A.x +B.x - C.x;
        C.fit_fun(ff,ud1, ud2);
        D.fit_fun(ff, ud1, ud2);
        for (int i = 0; i <= n - 3; i++)
        {
            if(C.y<D.y)
            {
                B=D;
            }
            else
            {
                A=C;
            }
            C.x = B.x - 1.0 * fib[n - i - 2] / fib[n -i -1] * (B.x - A.x);
            D.x = A.x+B.x - C.x;
            C.fit_fun(ff, ud1, ud2);
            D.fit_fun(ff, ud1, ud2);
        }
        return C;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
        static uint32_t call_count = 0;
        call_count++;
		solution Xopt;
        double A,B,C,D, prev_D;
		int i = 0;
        C =(A+B)/2;//srodek przedzialu
        do {
            matrix l = ff(A,0,0)*(pow(B,2)) - pow(C,2) - pow(A,2)+ ff(C,0,0)*(pow(a,2))- pow(b,2);
            matrix m = ff(A,0,0)*(B-C) + ff(B,0,0)*(C-A) + ff(C,0,0)*(A-B);
            if(m <= 0) throw;
            D = 0.5*det(l)/det(m);
            if(A<D<C)
            {
                if(ff(D,0,0)<ff(C,0,0))
                {
                    C=D;
                    B=C;
                }
                else
                {
                    A = D;
                }
            }
            else if (C<D<B)
            {
                if(ff(D,0,0)<ff(C,0,0))
                {
                    A = C;
                    C = D;
                }
                else
                {
                    B = D;
                }

            }
            else throw;
            i++;
            D = prev_D;
            if(call_count>Nmax) throw;
        }while((B-A)<epsilon || abs(D - prev_D)<gamma);

        Xopt = solution(D);
		return Xopt;
	}
	catch (string ex_info)
	{
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
                }
                while (X.y >= XB.y);
            }
            else
            {
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
