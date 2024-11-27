#include "Problem.h"

//переход от консервативных переменных к неконсервативным
void cons_to_noncons(Parameters& params, double& p, double& v_x, double& rho, double& m, double& imp_x, double& e) {
    p = (params.gamma - 1.0) * (e - 0.5 * pow(imp_x, 2.0) / m);
    v_x = imp_x / m;
    rho = m;
}

//переход от неконсервативных переменных к консервативным
void noncons_to_cons(Parameters& params, double& p, double& v_x, double& rho, double& m, double& imp_x, double& e) {
    m = rho;
    imp_x = rho * v_x;;
    e = 0.5 * rho * pow(v_x, 2.0) + p / (params.gamma - 1.0);
}

void Problem::init_grid() {
    for (int i = 0; i < params.N_x; i++) {
        x.push_back(params.start_x + (i - 1) * params.dx);
        xc.push_back(x[i] + params.dx / 2.);
    }
    x.push_back(params.start_x + (params.N_x - 1) * params.dx);
}

void Problem::init() {
    switch (params.initial) {
    case 1:
        Test1();
        break;
    case 2:
        Test2();
        break;
    case 3:
        Test3();
        break;
    case 4:
        Test4();
        break;
    case 5:
        Test5();
        break;
    default:
        Test1();
        break;
    }
    for (int i = 0; i < params.N_x - 1; i++) {
        ei[i] = p[i] / (rho[i] * (params.gamma - 1.));
        e[i] = ei[i] + pow(v_x[i], 2) / 2.;
        m[i] = (x[i + 1] - x[i]) * rho[i];
    }
}

void Problem::calc_xc() {
    for (int i = 0; i < params.N_x; i++) {
        xc[i] = (x[i] + x[i + 1]) / 2.;
    }
}

void Problem::get_dt() {
    double dt_min = 100;
    for (int i = 1; i < params.N_x - 1; ++i) {
        params.c = sqrt(params.gamma * p[i] / rho[i]);
        params.dt = (x[i + 1] - x[i]) / (abs(v_x[i]) + params.c);
        if (params.dt < dt_min) {
            dt_min = params.dt;
        }
    }
    params.dt = dt_min * params.CFL;
}

Problem::Problem()
{
    params = Parameters();
    init_grid();
}

Problem::Problem(Parameters& _params) : params(_params)
{
    init_grid();

    rho = vec(params.N_x);
    v_x = vec(params.N_x + 1);
    p = vec(params.N_x);
    pb = vec(params.N_x);

    rho_new = vec(params.N_x);
    v_x_new = vec(params.N_x + 1);
    p_new = vec(params.N_x);
    pb_new = vec(params.N_x);

    m = vec(params.N_x);
    imp_x = vec(params.N_x + 1);
    e = vec(params.N_x);
    ei = vec(params.N_x);

    x_new = vec(params.N_x + 1);
    imp_x_new = vec(params.N_x + 1);
    e_new = vec(params.N_x);
    ei_new = vec(params.N_x);

    w = vec(params.N_x);
}

Problem::~Problem()
{
}

void Problem::write_out(int _n, double _time)
{
    if (!fs::exists("result")) {
        fs::create_directory("result");
    }
    if (!fs::exists("result/x")) {
        fs::create_directory("result/x");
    }
    std::string file_path = "result/x/" + std::to_string(_n) + ".txt";
    std::ofstream csv;
    csv.open(file_path);
    csv << _time << std::endl;
    csv << "x;v_x;imp_x\n";
    for (int j = 1; j < params.N_x; j++) {
        csv << x[j] << ";" << v_x[j] << ";" << imp_x[j] << "\n";
    }
    csv.close();

    if (!fs::exists("result/xc")) {
        fs::create_directory("result/xc");
    }
    file_path = "result/xc/" + std::to_string(_n) + ".txt";
    csv.open(file_path);
    csv << _time << std::endl;
    csv << "xc;rho;m;p;e;ei\n";
    for (int j = 1; j < params.N_x - 1; j++) {
        csv << xc[j] << ";" << rho[j] << ";" << m[j] << ";" << p[j] << ";" << e[j] << ';' << ei[j] << "\n";
    }
    csv.close();
}

void Problem::Test1()
{
    double p1, rho1, v1;
    double p2, rho2, v2;
    rho1 = 1.0;
    v1 = 0.0;
    p1 = 1.0;
    rho2 = 0.125;
    v2 = 0.0;
    p2 = 0.1;
    for (int i = 0; i < params.N_x; ++i) {
        if (xc[i] < 0.5) {
            p[i] = p1;
            rho[i] = rho1;
        }
        else {
            p[i] = p2;
            rho[i] = rho2;
        }
    }
    for (int i = 0; i <= params.N_x; i++)
    {
        if (x[i] < 0.5)
        {
            v_x[i] = v1;
        }
        else
        {
            v_x[i] = v2;
        }
    }
}

void Problem::Test2()
{
    double p1, rho1, v1;
    double p2, rho2, v2;
    rho1 = 1.0;
    v1 = -2.0;
    p1 = 0.4;
    rho2 = 1.0;
    v2 = 2.0;
    p2 = 0.4;
    for (int i = 0; i < params.N_x; ++i) {
        if (xc[i] < 0.5) {
            p[i] = p1;
            rho[i] = rho1;
        }
        else {
            p[i] = p2;
            rho[i] = rho2;

        }
    }
    for (int i = 0; i <= params.N_x; i++)
    {
        if (x[i] < 0.5)
        {
            v_x[i] = v1;
        }
        else
        {
            v_x[i] = v2;
        }
    }
}

void Problem::Test3()
{
    double p1, rho1, v1;
    double p2, rho2, v2;
    rho1 = 1.0;
    v1 = 0.0;
    p1 = 1000.0;
    rho2 = 1.0;
    v2 = 0.0;
    p2 = 0.01;
    for (int i = 0; i < params.N_x; ++i) {
        if (xc[i] < 0.5) {
            p[i] = p1;
            rho[i] = rho1;
        }
        else {
            p[i] = p2;
            rho[i] = rho2;
        }
    }
    for (int i = 0; i <= params.N_x; i++)
    {
        if (x[i] < 0.5)
        {
            v_x[i] = v1;
        }
        else
        {
            v_x[i] = v2;
        }
    }
}

void Problem::Test4() {
    double p1, rho1, v1;
    double p2, rho2, v2;
    rho1 = 1.0;
    v1 = 0.0;
    p1 = 0.01;
    rho2 = 1.0;
    v2 = 0.0;
    p2 = 100.0;
    for (int i = 0; i < params.N_x; ++i) {
        if (xc[i] < 0.5) {
            p[i] = p1;
            rho[i] = rho1;
        }
        else {
            p[i] = p2;
            rho[i] = rho2;
        }
    }
    for (int i = 0; i <= params.N_x; i++)
    {
        if (x[i] < 0.5)
        {
            v_x[i] = v1;
        }
        else
        {
            v_x[i] = v2;
        }
    }
}

void Problem::Test5() {
    double p1, rho1, v1;
    double p2, rho2, v2;
    rho1 = 5.99924;
    v1 = 19.5975;
    p1 = 460.894;
    rho2 = 5.99242;
    v2 = -6.19633;
    p2 = 46.0950;
    for (int i = 0; i < params.N_x; ++i) {
        if (xc[i] < 0.5) {
            p[i] = p1;
            rho[i] = rho1;
        }
        else {
            p[i] = p2;
            rho[i] = rho2;
        }
    }
    for (int i = 0; i <= params.N_x; i++)
    {
        if (x[i] < 0.5)
        {
            v_x[i] = v1;
        }
        else
        {
            v_x[i] = v2;
        }
    }
}

void Problem::ArtificialViscosity()
{
    double mu_0 = 2.;
    for (int i = 1; i < params.N_x - 1; ++i) {
        if (v_x[i + 1] - v_x[i] < 0.) {
            w[i] = -mu_0 * rho[i] * abs(v_x[i + 1] - v_x[i]) * (v_x[i + 1] - v_x[i]);
        }
        else {
            w[i] = 0.;
        }
    }
    w[0] = w[1];
    w[params.N_x - 1] = w[params.N_x - 2];
    for (int i = 0; i < params.N_x; ++i) {
        p[i] = p[i] + w[i];
    }
}


void Problem::ArtificialViscosityConst()
{
    double mu_0 = 2.;
    for (int i = 1; i < params.N_x - 1; ++i) {
        if (v_x[i + 1] - v_x[i] < 0.) {
            w[i] = -mu_0 * m[i];
        }
        else {
            w[i] = 0.;
        }
    }
    w[0] = w[1];
    w[params.N_x - 1] = w[params.N_x - 2];
    for (int i = 0; i < params.N_x; ++i) {
        p[i] = p[i] + w[i];
    }
}

void Problem::ArtificialViscosityLinear()
{
    double mu_0 = 2.;
    for (int i = 1; i < params.N_x - 1; ++i) {
        w[i] = -mu_0 * rho[i] * (v_x[i + 1] - v_x[i]);
    }
    w[0] = w[1];
    w[params.N_x - 1] = w[params.N_x - 2];
    for (int i = 0; i < params.N_x; ++i) {
        p[i] = p[i] + w[i];
    }
}

void Problem::ArtificialViscosityNeumann()
{
    double mu_0 = 2.;
    for (int i = 1; i < params.N_x - 1; ++i) {
        w[i] = -mu_0 * rho[i] * rho[i] * abs(v_x[i + 1] - v_x[i]) * (v_x[i + 1] - v_x[i]);
    }
    w[0] = w[1];
    w[params.N_x - 1] = w[params.N_x - 2];
    for (int i = 0; i < params.N_x; ++i) {
        p[i] = p[i] + w[i];
    }
}

void Problem::ArtificialViscosityQuadro()
{
    double mu_0 = 2.;
    for (int i = 1; i < params.N_x - 1; ++i) {
        if (v_x[i + 1] - v_x[i] < 0.) {
            w[i] = mu_0 * rho[i] * (v_x[i + 1] - v_x[i]) * (v_x[i + 1] - v_x[i]);
        }
        else {
            w[i] = 0.;
        }
    }
    w[0] = w[1];
    w[params.N_x - 1] = w[params.N_x - 2];
    for (int i = 0; i < params.N_x; ++i) {
        p[i] = p[i] + w[i];
    }
}

void Problem::ArtificialViscosityMixed()
{
    double mu_0 = 2.;
    for (int i = 1; i < params.N_x - 1; ++i) {
        if (v_x[i + 1] - v_x[i] < 0.) {
            w[i] = mu_0 * rho[i] * (v_x[i + 1] - v_x[i]) * (v_x[i + 1] - v_x[i]) - mu_0 * rho[i] * (v_x[i + 1] - v_x[i]);
        }
        else {
            w[i] = 0.;
        }
    }
    w[0] = w[1];
    w[params.N_x - 1] = w[params.N_x - 2];
    for (int i = 0; i < params.N_x; ++i) {
        p[i] = p[i] + w[i];
    }
}

void Problem::Solve()
{
    double time = params.start_t;
    int n = 0;
    init();
    if (fs::exists("result")) {
        fs::remove_all("result");
    }
    write_out(n, time);
    while (time < params.end_t && n < params.max_iter)
    {
        get_dt();
        n++;
        time += params.dt;
        ArtificialViscosityMixed();
        if (params.solver == 2) {
            CrossSolverDivergent();
        }
        else {
            CrossSolverNonDivergent();
        }
        if (n % params.write_interval == 0) {
            write_out(n, time);
        }
    }
}

void Problem::CrossSolverNonDivergent()
{
    for (int j = 1; j < params.N_x; j++)
    {
        v_x_new[j] = v_x[j] - params.dt / (0.5 * (m[j] + m[j - 1])) * (p[j] - p[j - 1]);
    }

    for (int j = 1; j < params.N_x; j++)
    {
        x_new[j] = x[j] + params.dt * v_x_new[j];
        imp_x_new[j] = v_x[j] * m[j];
    }

    for (int j = 1; j < params.N_x - 1; j++)
    {
        rho_new[j] = 1. / (params.dt / m[j] * (v_x_new[j + 1] - v_x_new[j]) + 1. / rho[j]);
        //rho_new[j] = m[j] / (x_new[j + 1] - x_new[j]);
        ei_new[j] = ei[j] / (1. + rho_new[j] * (params.gamma - 1.) * params.dt / m[j] * (v_x_new[j + 1] - v_x_new[j]));
        e_new[j] = ei_new[j] + pow((v_x[j] + v_x[j + 1] / 2.), 2) / 2.;
        p_new[j] = rho_new[j] * (params.gamma - 1.) * ei_new[j];
    }

    v_x_new[1] = 0.;
    v_x_new[params.N_x - 1] = 0.;
    p_new[0] = p_new[1];
    rho_new[0] = rho_new[1];
    p_new[params.N_x - 1] = p_new[params.N_x - 2];
    rho_new[params.N_x - 1] = rho_new[params.N_x - 2];

    x = x_new;
    calc_xc();
    v_x = v_x_new;
    rho = rho_new;
    p = p_new;
    e = e_new;
    ei = ei_new;
}

void Problem::CrossSolverDivergent()
{

    for (int j = 1; j < params.N_x; j++)
    {
        pb[j] = (p[j] * m[j - 1] + p[j - 1] * m[j]) / (m[j] + m[j - 1]);
    }
    pb[0] = pb[1];
    pb[params.N_x - 1] = pb[params.N_x - 2];

    for (int j = 1; j < params.N_x; j++)
    {
        v_x_new[j] = v_x[j] - params.dt * (p[j] - p[j - 1]) / (0.5 * (m[j] + m[j - 1]));
    }

    for (int j = 1; j < params.N_x; j++)
    {
        x_new[j] = x[j] + params.dt * v_x_new[j];
        imp_x_new[j] = v_x[j] * m[j];
    }

    for (int j = 1; j < params.N_x - 1; j++)
    {
        rho_new[j] = m[j] / (x_new[j + 1] - x_new[j]);
        e_new[j] = e[j] - params.dt * (pb[j + 1] * v_x_new[j + 1] - pb[j] * v_x_new[j]) / m[j];
        ei_new[j] = e_new[j] - pow(((v_x[j] + v_x[j + 1]) / 2.), 2) / 2.;
        p_new[j] = (params.gamma - 1) * ei_new[j] * rho_new[j];
    }

    v_x_new[1] = 0.;
    v_x_new[params.N_x - 1] = 0.;
    p_new[0] = p_new[1];
    rho_new[0] = rho_new[1];
    p_new[params.N_x - 1] = p_new[params.N_x - 2];
    rho_new[params.N_x - 1] = rho_new[params.N_x - 2];

    x = x_new;
    calc_xc();
    v_x = v_x_new;
    rho = rho_new;
    p = p_new;
    e = e_new;
    ei = ei_new;
}