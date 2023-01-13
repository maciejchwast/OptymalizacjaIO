#include"opt_alg.h"

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x, double d, double alpha, int N_max, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0,0 };
        int i = 0;
        solution X0(x), X1(x + d);
        X0.fit_fun(ff,ud1,ud2);
        X1.fit_fun(ff,ud1,ud2);
        if(X0.y == X1.y){
            p[0] =m2d(X0.x);
            p[1] = m2d(X1.x);
            return p;
        }

        if(X0.y < X1.y){
            d = -d;
            X1.x = X0.x + d;
            X1.fit_fun(ff,ud1,ud2);
            if(X1.y >= X0.y){
                p[0] = m2d(X1.x);
                p[1] = m2d(X0.x)-d;
                return p;
            }
        }
        solution X2;
        while(true)
        {
            i++;
            X2.x = x + pow(alpha,i)*d;
            X2.fit_fun(ff,ud1,ud2);
            if(solution::f_calls>N_max || X2.y >= X1.y){
                break;
            }
            X0=X1;
            X1=X2;
        }

        if(d>0)
        {
            p[0] = m2d(X0.x);
            p[1] = m2d(X2.x);
            return p;
        }
        p[0] = m2d(X2.x);
        p[1] = m2d(X0.x);
        return p;

	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		Xopt.ud = b-a;
        int n = static_cast<int>(ceil(log2(sqrt(5) * (b - a) / epsilon) / log2((1 + sqrt(5)) / 2)));
        int* F = new int[n] {1, 1};
        for (int i = 2; i < n; ++i)
            F[i] = F[i - 2] + F[i - 1];
        solution A(a), B(b), C, D;
        C.x = B.x - 1.0 * F[n - 2] / F[n - 1] * (B.x - A.x);
        D.x = A.x + B.x - C.x;
        C.fit_fun(ff, ud1, ud2);
        D.fit_fun(ff, ud1, ud2);
        for (int i = 0; i <= n - 3; ++i) {
            if(C.y< D.y)
            {
                B=D;
            }
            else
            {
                A = C;
            }
            C.x = B.x - 1.0 * F[n-i-2] / F[n-i-1] * (B.x-A.x);
            D.x = A.x + B.x - C.x;
            C.fit_fun(ff,ud1,ud2);
            D.fit_fun(ff,ud1,ud2);

            Xopt.ud.add_row((B.x-A.x)());
        }

        Xopt = C;
        Xopt.flag = 0;
        return Xopt;

        /*for(int i=2;i<100;i++)
        {
            fibNumbers.push_back(fibNumbers.at(i-1) + fibNumbers.at(i-2));
        }

        int k = 100;
        for(int i=0;i<fibNumbers.size();i++)
        {
            if(fibNumbers.at(i) > (b - a) / epsilon){
                k = i;
                break;
            }
        }

        solution A, B, C, D;

        A = a;
        B = b;
        C.x = B.x - fibNumbers.at(k-1) / fibNumbers.at(k) * (B.x-A.x);
        D.x = A.x + B.x - C.x;
        C.fit_fun(ff,ud1,ud2);
        D.fit_fun(ff,ud1,ud2);
        for(int i=0;i<=k-3;i++)
        {
            if(fun1(C.y)< fun1(D.y))
            {
                B=D;
            }
            else
            {
                A = C;
            }
            C.x = B.x - fibNumbers.at(k - i - 2) / fibNumbers.at(k-i-1) * (B.x-A.x);
            D.x = A.x + B.x - C.x;
            C.fit_fun(ff,ud1,ud2);
            D.fit_fun(ff,ud1,ud2);
        }

        Xopt = solution(C);
		return Xopt;*/
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
        solution Xopt;
        Xopt.ud = b - a;

        solution A(a), B(b), C, D, prev_D(a);
        C.x = (a + b) / 2;
        A.fit_fun(ff, ud1, ud2);
        B.fit_fun(ff, ud1, ud2);
        C.fit_fun(ff, ud1, ud2);
        double l, m;
        while(true)
        {
            l = m2d(A.y * (pow(B.x) - pow(C.x)) + B.y * (pow(C.x) - pow(A.x)) + C.y *(pow(A.x) - pow(B.x)));
            m = m2d(A.y * (B.x - C.x) + B.y * (C.x - A.x) + C.y * (A.x - B.x));
            if (m <= 0)
            {
                Xopt = prev_D;
                Xopt.flag = 2;
                return Xopt;
            }
            D.x = 0.5 * l / m;
            D.fit_fun(ff, ud1, ud2);

            if (A.x <= D.x && D.x <= C.x)
            {
                if(D.y<C.y)
                {
                    B=C;
                    C=D;
                }
                else
                {
                   A=D;
                }

            }
            else if (C.x <= D.x && D.x <= B.x)
            {
                if(D.y<C.y)
                {
                    A.x = C.x;
                    C.x = D.x;
                }
                else
                {
                    B=D;

                }
            }
            else
            {
                Xopt=prev_D;
                Xopt.flag=2;
                return Xopt;
            }

            Xopt.ud.add_row((B.x-A.x)());

            if(B.x-A.x<epsilon || abs(D.x() - prev_D.x())<gamma)
            {
                A.fit_fun(ff,ud1,ud2);
                B.fit_fun(ff,ud1,ud2);
                C.fit_fun(ff,ud1,ud2);

                Xopt=D;
                Xopt.flag=0;
                break;
            }
            if(solution::f_calls>Nmax)
            {
                A.fit_fun(ff,ud1,ud2);
                B.fit_fun(ff,ud1,ud2);
                C.fit_fun(ff,ud1,ud2);

                Xopt = D;
                Xopt.flag = 1;
                break;
            }
            A.fit_fun(ff,ud1,ud2);
            B.fit_fun(ff,ud1,ud2);
            C.fit_fun(ff,ud1,ud2);
            prev_D=D;
        }
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
        solution XB(x0), XB_old, X;
        XB.fit_fun(ff,ud1, ud2);
        while (true)
        {
            X = HJ_trial(ff,XB, s, ud1, ud2);
            if (X.y<XB.y) //etapy robocze wywolujemy tak dlugo az przynosza poprawe
            {
                while (true)
                {
                    XB_old = XB;
                    XB = X;
                    X.x = 2*XB.x -XB_old.x;
                    X.fit_fun(ff,ud1, ud2);
                    X = HJ_trial(ff, X, s, ud1, ud2);
                    if (X.y>=XB.y) //jezeli jest gorszy od x bazowe
                        break;
                    if (solution::f_calls>Nmax) //czy nie przekroczylismy liczby iteracji
                        return XB;
                }
            }
            else
                s*=alpha; //je�eli nie ma poprawy zmniejszamy krok
            if (s<epsilon || solution::f_calls>Nmax)
                return XB;
        }
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
        int n = get_dim(XB);
        matrix E = ident_mat(n);
        solution X;
        for(int j = 0; j < n; j++){
            X.x = XB.x+s*E[j];
            X.fit_fun(ff,ud1,ud2);
            if(X.y < XB.y)
            {
                XB = X;
            }
            else
            {
                X.x = XB.x - s*E[j];
                X.fit_fun(ff,ud1,ud2);
                if(X.y < XB.y)
                {
                    XB = X;
                }
            }
        }
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
        solution X(x0), Xt;
        int n = get_dim(X);
        matrix l(n, 1), p(n, 1), s(s0), D = ident_mat(n); // l - lambda, p - licznik porazek, s - dlugosci krokow, d - macierz kierunkami poczatkowtmi
        X.fit_fun(ff,ud1,ud2);
        while (true)
        {
            for (int i = 0; i < n; i++)
            {
                Xt.x = X.x + s(i) * D[i];
                Xt.fit_fun(ff,ud1,ud2);
                if (Xt.y <X.y)
                {
                    X = Xt;
                    l(i) += s(i);
                    s(i) *= alpha;
                }
                else
                {
                    s(i) *= -beta;
                    ++p(i);
                }
            }

            bool change = true;
            for (int i = 0; i < n; ++i)
                if (p(i)==0 || l(i)==0)
                {
                    change = false;
                    break;
                }
            if (change)
            {
                matrix Q(n,n), v(n,1);
                for (int i = 0; i<n; ++i)
                    for (int j = 0; j<=i; ++j)
                        Q(i, j) = l(i);
                Q = D*Q;
                v = Q[0] / norm(Q[0]);
                D.set_col(v,0);
                for (int i = 1; i < n; ++i)
                {
                    matrix temp(n,1);
                    for (int j = 0; j < i; ++j)
                        temp = temp + trans(Q[i]) * D[j]*D[j];

                    v = (Q[i] - temp) / norm(Q[i] - temp);
                    D.set_col(v,i);
                }
                s = s0;
                p = matrix(n, 1);
                l = matrix(n, 1);
            }
            double max_s = abs(s(0));
            for (int i = 1; i < n; ++i)
                if (max_s < abs(s(i)))
                    max_s = abs(s(i));
            if (max_s < epsilon || solution::f_calls > Nmax)
                return X;
        }
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c0, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
        double alpha = 1, beta = 0.5, gamma = 2, delta = 0.5, s = 0.5;
        solution X(x0), X1;
        matrix c(2, new double[2]{c0, dc});
        while (true) {
            X1 = sym_NM(ff,X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c);
            if (norm(X1.x - X.x) < epsilon || solution::f_calls > Nmax)
                return X1;
            c(0) = dc * c(0); //jezeli nie konczymy zmieniamy wspolczynnik c
            X = X1;
        }
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
        int n = get_len(x0);
        matrix D = ident_mat(n);
        int N = n + 1;
        solution *S = new solution[N];
        S[0].x = x0;
        S[0].fit_fun(ff,ud1, ud2);
        for (int i = 1; i < N; ++i)
        {
            S[i].x = S[0].x + s*D[i-1];
            S[i].fit_fun(ff,ud1, ud2);
        }
        solution PR, PE, PN; //PR - odbite, PE - ekspansja, PN - zawezenie
        matrix pc; //srodek ciezkosci simpleksu
        int i_min, i_max; //najlepszy i najgorszy wiercholek
        while (true)
        {
            i_min = i_max = 0;
            for (int i = 1; i < N; ++i)
            {
                if (S[i_min].y>S[i].y)
                    i_min = i;
                if (S[i_max].y<S[i].y)
                    i_max = i;
            }
            pc = matrix(n,1);
            for (int i = 0; i < N; ++i)
                if (i!=i_max)
                    pc=pc+S[i].x;
            pc = pc / (N-1);
            PR.x = pc + alpha*(pc - S[i_max].x);
            PR.fit_fun(ff,ud1, ud2);
            if (PR.y < S[i_max].y && S[i_min].y <= PR.y)
                S[i_max] = PR;
            else if (PR.y<S[i_min].y)
            {
                PE.x = pc+gamma*(PR.x-pc);
                PE.fit_fun(ff,ud1, ud2);
                if (PE.y < S[i_max].y)
                    S[i_max] = PE;
                else
                    S[i_max] = PR;
            }
            else
            {
                PN.x = pc+beta*(S[i_max].x-pc);
                PN.fit_fun(ff,ud1, ud2);
                if (PN.y < S[i_max].y)
                    S[i_max] = PN;
                else
                {
                    for (int i = 0; i < N; ++i)
                        if (i!=i_min)
                        {
                            S[i].x = delta*(S[i].x+S[i_min].x);
                            S[i].fit_fun(ff,ud1, ud2);
                        }
                }
            }

            double max_s = norm(S[0].x - S[i_min].x);
            for (int i = 1; i < N; ++i)
                if (norm(S[i].x-S[i_min].x)>max_s)
                    max_s = norm(S[i].x - S[i_min].x);
            if (max_s < epsilon || solution::f_calls>Nmax)
                return S[i_min];
        }
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
        int n = get_len(x0);
        solution X, X1;
        X.x = x0;
        matrix d(n, 1), *P = new matrix[2]; //macierz kierunku d, p: dwuelementowa macierz
        solution h;
        double *ab;
        while (true) {
            X.grad(gf);
            d = -X.g;
            if (h0 < 0) {
                P[0] = X.x;
                P[1] = d;
                ab = expansion(ff,0, 1, 1.2, Nmax, ud1, *P);
                h = golden(ff,ab[0], ab[1], epsilon, Nmax, ud1, *P);
                X1.x = X.x - h.x * d;
            } else
                X1.x = X.x + h0 * d;



            if (norm(X1.x - X.x) < epsilon || solution::f_calls > Nmax || solution::g_calls > Nmax) {
                X1.fit_fun(ff,ud1, ud2);
                return X1;
            }
            X = X1;
        }
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
        int n = get_len(x0);
        solution X, X1;
        X.x = x0;
        matrix d(n, 1), *P = new matrix[2];
        solution h;
        double *ab, beta;
        X.grad(gf);
        d = -X.g;
        while (true)
        {
            if (h0<0)
            {
                P[0] = X.x;
                P[1] = d;
                ab = expansion(fun4,0,0.2,1.2,Nmax, ud1, *P);
                h = golden(fun4,ab[0], ab[1], epsilon, Nmax, ud1,*P);
                X1.x = X.x + h.x * d;
            }
            else
                X1.x = X.x + h0 * d;

            if (norm(X1.x - X.x)<epsilon || solution::f_calls>Nmax || solution::g_calls>Nmax)
            {
                X1.fit_fun(ff,ud1);
                return X1;
            }
            X1.grad(gf);
            beta = pow(norm(X1.g),2)/ pow(norm(X.g),2);
            d = -X1.g+beta*d;
            X =X1;
        }
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),matrix(*hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
        int n = get_len(x0);
        solution X, X1;
        X.x = x0;
        matrix d(n, 1), *P = new matrix[2];
        solution h;
        double *ab;
        while (true)
        {
            X.grad(gf);
            X.hess(hf);
            d = -inv(X.H)*X.g; // inv macierz odwrotna
            if (h0<0)
            {
                P[0] = X.x;
                P[1] = d;
                ab = expansion(ff,0,1,1.2,Nmax, ud1, *P);
                h = golden(ff,ab[0], ab[1], epsilon, Nmax, ud1,*P);
                X1.x = X.x + h.x * d;
            }
            else
                X1.x = X.x + h0 * d;


            if (norm(X1.x - X.x)<epsilon || solution::f_calls>Nmax || solution::g_calls>Nmax)
            {
                X1.fit_fun(ff,ud1);
                return X1;
            }
            X=X1;
        }
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
        double alfa = (sqrt(5)-1)/2;
        solution A, B, C, D;
        A.x = a;
        B.x = b;
        C.x = B.x - alfa*(B.x-A.x);
        C.fit_fun(ff,ud1, ud2);
        D.x = A.x + alfa*(B.x-A.x);
        D.fit_fun(ff,ud1, ud2);
        while (true)
        {
            if (C.y<D.y)
            {
                B=D;
                D=C;
                C.x = B.x - alfa*(B.x-A.x);
                C.fit_fun(ff,ud1, ud2);
            }
            else
            {
                A = C;
                C=D;
                D.x = A.x + alfa*(B.x-A.x);
                D.fit_fun(ff,ud1, ud2);
            }
            if (B.x-A.x<epsilon || solution::f_calls > Nmax)
            {
                A.x = (A.x+B.x)/2;
                A.fit_fun(ff,ud1, ud2);
                return A;
            }
        }
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
		int n = get_len(x0); //dlugosc wektora
        matrix D = ident_mat(n), *A = new matrix[2]; //macierz kierunkow, m. jednostkowa
        //A-wskaznik, A0 -punkt gdzie jestesmy, A1-kierunek
        solution X, P, h; //X-punkt, p-p0,p1...pn,h-dlugosc kroku
        X.x=x0;
        double *ab;
        while(true)
        {
            P=X; //start od punktu x
            for(int i =0; i<n; i++)
            {
                    A[0] = P.x; //punkt
                    A[1] = D[i]; //kierunek
                    ab = expansion(ff,0,1,1.2,Nmax,ud1,*A);
                    h = golden(ff,ab[0],ab[1],epsilon,Nmax,ud1,*A);
                    P.x = P.x +h.x *D[i]; //przesuniecie
            }
            if(norm(P.x - X.x) <epsilon || solution::f_calls>Nmax)
            {
                P.fit_fun(ff,ud1,ud2);
                return P;
            }
            for(int i = 0; i<n-1; i++)
            {
                D.set_col(D[i+1],i); //w miejsce kolumny i wstawiamy kolumne i+1
            }
            D.set_col(P.x - X.x, n-1);
            A[0] = P.x;
            A[1] = D[n-1];
            ab = expansion(ff, 0,1,1.2,Nmax, ud1,*A);
            h = golden(ff,ab[0],ab[1],epsilon, Nmax,ud1,*A);
            X.x = P.x + h.x *D[n-1]; //zmiana punktu startowego
        }
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
