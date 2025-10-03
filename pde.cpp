#include <bits/stdc++.h>
using namespace std;

/* ----------------- Types & Params ----------------- */
enum class PDEType { HEAT, WAVE, LAPLACE };
enum class Scheme { EXPLICIT, IMPLICIT, CRANK_NICOLSON };

struct PDEParams {
    double L = 1.0;
    int Nx = 101;
    double T = 0.1;
    double alpha = 0.01;   // diffusion coeff for heat
    double c = 1.0;        // wave speed
    double dt = 0.0;       // if 0 -> auto choose
    int Nt = 0;            // if >0 overrides dt by Nt
};

void write_csv(const string &fname, const vector<double>& x, const vector<double>& u) {
    ofstream f(fname);
    f << "x,u\n";
    f << fixed << setprecision(12);
    for (size_t i = 0; i < x.size(); ++i) f << x[i] << "," << u[i] << "\n";
    f.close();
}

vector<double> thomas_solve(const vector<double>& a, const vector<double>& b,
                            const vector<double>& c, const vector<double>& d) {
    // Solve tridiagonal system: a_i * x_{i-1} + b_i * x_i + c_i * x_{i+1} = d_i
    // a[0] is unused (or zero), c[n-1] unused (or zero)
    int n = (int)b.size();
    vector<double> cp(n, 0.0), dp(n, 0.0), x(n, 0.0);
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for (int i = 1; i < n; ++i) {
        double den = b[i] - a[i] * cp[i-1];
        cp[i] = (i < n-1) ? c[i] / den : 0.0;
        dp[i] = (d[i] - a[i] * dp[i-1]) / den;
    }
    x[n-1] = dp[n-1];
    for (int i = n-2; i >= 0; --i) x[i] = dp[i] - cp[i] * x[i+1];
    return x;
}

/* ----------------- Heat: (FTCS) ----------------- */
void heat_explicit(const PDEParams& p,
    function<double(double)> ic,
    function<double(double)> bcL,
    function<double(double)> bcR,
    const string &outfile = "heat_explicit.csv")
{
    int Nx = p.Nx;
    double L = p.L;
    double dx = L / (Nx - 1);
    double alpha = p.alpha;

    // choose dt through stability if not provided
    double dt;
    if (p.Nt > 0) { dt = p.T / p.Nt; }
    else if (p.dt > 0) { dt = p.dt; }
    else { dt = 0.4 * dx*dx / alpha; } // stable: dt <= 0.5*dx^2/alpha

    int Nt = max(1, (int)round(p.T / dt));
    dt = p.T / Nt;

    vector<double> u(Nx), u_new(Nx);
    vector<double> x(Nx);
    for (int i = 0; i < Nx; ++i) { x[i] = i * dx; u[i] = ic(x[i]); }

    for (int n = 0; n < Nt; ++n) {
        double t = n * dt;
        double tnew = (n+1) * dt;
        // interior updates
        double r = alpha * dt / (dx*dx);
        for (int i = 1; i < Nx-1; ++i)
            u_new[i] = u[i] + r * (u[i+1] - 2*u[i] + u[i-1]);
        // Dirichlet BCs at new time
        u_new[0] = bcL(tnew);
        u_new[Nx-1] = bcR(tnew);
        u.swap(u_new);
    }
    write_csv(outfile, x, u);
    cout << "heat_explicit -> " << outfile << "\n";
}

/* ----------------- Heat: (Backward Euler) ----------------- */
void heat_implicit(const PDEParams& p,
    function<double(double)> ic,
    function<double(double)> bcL,
    function<double(double)> bcR,
    const string &outfile = "heat_implicit.csv")
{
    int Nx = p.Nx;
    double L = p.L;
    double dx = L / (Nx - 1);
    double alpha = p.alpha;

    double dt;
    if (p.Nt > 0) { dt = p.T / p.Nt; }
    else if (p.dt > 0) { dt = p.dt; }
    else { dt = p.T / 200.0; } // default time resolution for implicit
    int Nt = max(1, (int)round(p.T / dt));
    dt = p.T / Nt;

    vector<double> u(Nx), u_new(Nx), x(Nx);
    for (int i = 0; i < Nx; ++i) { x[i] = i*dx; u[i] = ic(x[i]); }

    int ninternal = Nx - 2;
    vector<double> a(ninternal, 0.0), b(ninternal, 0.0), c(ninternal, 0.0), rhs(ninternal, 0.0);

    double r = alpha * dt / (dx*dx);
    // coefficients (same each time)
    for (int i = 0; i < ninternal; ++i) {
        a[i] = -r; b[i] = 1 + 2*r; c[i] = -r;
    }
    a[0] = 0.0; c[ninternal-1] = 0.0;

    for (int n = 0; n < Nt; ++n) {
        double tnew = (n+1) * dt;
        // RHS = u_old_interior
        for (int i = 0; i < ninternal; ++i) rhs[i] = u[i+1];
        // incorporate Dirichlet BCs at new time
        rhs[0] += r * bcL(tnew);
        rhs[ninternal-1] += r * bcR(tnew);

        // solve
        vector<double> sol = thomas_solve(a,b,c,rhs);
        // place solution
        u_new[0] = bcL(tnew);
        for (int i = 0; i < ninternal; ++i) u_new[i+1] = sol[i];
        u_new[Nx-1] = bcR(tnew);
        u.swap(u_new);
    }
    write_csv(outfile, x, u);
    cout << "heat_implicit -> " << outfile << "\n";
}

/* ----------------- Heat: Crank-Nicolson ----------------- */
void heat_crank_nicolson(const PDEParams& p,
    function<double(double)> ic,
    function<double(double)> bcL,
    function<double(double)> bcR,
    const string &outfile = "heat_cn.csv")
{
    int Nx = p.Nx;
    double L = p.L;
    double dx = L / (Nx - 1);
    double alpha = p.alpha;

    double dt;
    if (p.Nt > 0) { dt = p.T / p.Nt; }
    else if (p.dt > 0) { dt = p.dt; }
    else { dt = p.T / 200.0; }
    int Nt = max(1, (int)round(p.T / dt));
    dt = p.T / Nt;

    vector<double> u(Nx), u_new(Nx), x(Nx);
    for (int i = 0; i < Nx; ++i) { x[i] = i*dx; u[i] = ic(x[i]); }

    int ninternal = Nx - 2;
    vector<double> a(ninternal, 0.0), b(ninternal, 0.0), c(ninternal, 0.0), rhs(ninternal, 0.0);

    double r = alpha * dt / (dx*dx);

    for (int i = 0; i < ninternal; ++i) {
        a[i] = -r/2.0;
        b[i] = 1.0 + r;
        c[i] = -r/2.0;
    }
    a[0] = 0.0; c[ninternal-1] = 0.0;

    for (int n = 0; n < Nt; ++n) {
        double told = n * dt;
        double tnew = (n+1) * dt;
        // RHS = u_old + (r/2) * lap(u_old)
        for (int i = 0; i < ninternal; ++i) {
            int gi = i+1;
            rhs[i] = u[gi] + (r/2.0) * (u[gi+1] - 2.0*u[gi] + u[gi-1]);
        }

        rhs[0] += (r/2.0) * bcL(tnew);
        rhs[ninternal-1] += (r/2.0) * bcR(tnew);
        
        vector<double> sol = thomas_solve(a,b,c,rhs);
        u_new[0] = bcL(tnew);
        for (int i = 0; i < ninternal; ++i) u_new[i+1] = sol[i];
        u_new[Nx-1] = bcR(tnew);
        u.swap(u_new);
    }
    write_csv(outfile, x, u);
    cout << "heat_crank_nicolson -> " << outfile << "\n";
}

/* ----------------- Wave: Explicit ----------------- */
void wave_explicit(const PDEParams& p,
    function<double(double)> ic,
    function<double(double)> vel0,           
    function<double(double)> bcL,
    function<double(double)> bcR,
    const string &outfile = "wave_explicit.csv")
{
    int Nx = p.Nx;
    double L = p.L;
    double dx = L / (Nx - 1);
    double c = p.c;

    double dt;
    if (p.Nt > 0) { dt = p.T / p.Nt; }
    else if (p.dt > 0) { dt = p.dt; }
    else { dt = 0.9 * dx / c; } 

    int Nt = max(1, (int)round(p.T / dt));
    dt = p.T / Nt;

    vector<double> u_prev(Nx), u_curr(Nx), u_next(Nx), x(Nx);
    for (int i = 0; i < Nx; ++i) { x[i] = i*dx; u_prev[i] = ic(x[i]); }

    double t0 = 0.0;
    double t1 = dt;
    for (int i = 1; i < Nx-1; ++i) {
        double lap = (u_prev[i+1] - 2*u_prev[i] + u_prev[i-1])/(dx*dx);
        u_curr[i] = u_prev[i] + dt * vel0(x[i]) + 0.5 * (c*c) * (dt*dt) * lap;
    }
    u_curr[0] = bcL(t1);
    u_curr[Nx-1] = bcR(t1);

    double r = (c*c) * (dt*dt) / (dx*dx);

    for (int n = 1; n < Nt; ++n) {
        double t = n * dt;
        double tnew = (n+1) * dt;
        for (int i = 1; i < Nx-1; ++i) {
            u_next[i] = 2.0*u_curr[i] - u_prev[i] + r * (u_curr[i+1] - 2.0*u_curr[i] + u_curr[i-1]);
        }
        u_next[0] = bcL(tnew);
        u_next[Nx-1] = bcR(tnew);
        u_prev.swap(u_curr);
        u_curr.swap(u_next);
    }
    write_csv(outfile, x, u_curr);
    cout << "wave_explicit -> " << outfile << "\n";
}

/* ----------------- Laplace / Poisson: Jacobi ----------------- */
void laplace_jacobi(const PDEParams& p,
    function<double(double)> f,  
    function<double(double)> bcL,
    function<double(double)> bcR,
    const string &outfile = "laplace.csv",
    int maxIter = 10000,
    double tol = 1e-8)
{
    int Nx = p.Nx;
    double L = p.L;
    double dx = L / (Nx - 1);

    vector<double> u(Nx), u_new(Nx), x(Nx);
    for (int i = 0; i < Nx; ++i) { x[i] = i*dx; u[i] = 0.0; }
    u[0] = bcL(0.0); u[Nx-1] = bcR(0.0);
    for (int i = 1; i < Nx-1; ++i) {
        u[i] = ( (Nx-1-i) * u[0] + i * u[Nx-1] ) / (Nx-1);
    }

    for (int iter = 0; iter < maxIter; ++iter) {
        double maxdiff = 0.0;
        for (int i = 1; i < Nx-1; ++i) {
            double xi = x[i];
            double rhs = - (dx*dx) * f(xi); // scheme: u_{i-1} - 2u_i + u_{i+1} = dx^2 * f -> u_i = 0.5*(u_{i-1}+u_{i+1} - dx^2 f)
            u_new[i] = 0.5 * (u[i-1] + u[i+1] + rhs);
            maxdiff = max(maxdiff, fabs(u_new[i] - u[i]));
        }
        for (int i = 1; i < Nx-1; ++i) u[i] = u_new[i];
        if (maxdiff < tol) break;
    }
    write_csv(outfile, x, u);
    cout << "laplace_jacobi -> " << outfile << "\n";
}

int main() {
    cout << "PDE MVP demo\n";
    PDEParams p;
    p.L = 1.0;
    p.Nx = 101;
    p.T = 0.1;
    p.alpha = 0.01;

   
    auto ic = [&](double x){ return sin(M_PI * x); };
    auto zeroBC = [&](double t){ return 0.0; };

    heat_explicit(p, ic, zeroBC, zeroBC, "heat_explicit.csv");
    heat_implicit(p, ic, zeroBC, zeroBC, "heat_implicit.csv");
    heat_crank_nicolson(p, ic, zeroBC, zeroBC, "heat_cn.csv");

    PDEParams pw = p;
    pw.T = 0.5;
    pw.c = 1.0;
    auto vel0 = [&](double x){ return 0.0; };
    wave_explicit(pw, ic, vel0, zeroBC, zeroBC, "wave_explicit.csv");

    // Laplace: f(x)=0 (pure Laplace), BCs 0
    PDEParams pl = p;
    pl.T = 0.0;
    auto fzero = [&](double x){ return 0.0; };
    laplace_jacobi(pl, fzero, zeroBC, zeroBC, "laplace.csv");

    cout << "Demo complete. CSV files: heat_explicit.csv, heat_implicit.csv, heat_cn.csv, wave_explicit.csv, laplace.csv\n";
    return 0;
}
