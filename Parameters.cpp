#include "Parameters.h"

void Parameters::calc_dt()
{
    dt = 0.01;
}

void Parameters::calc_dx()
{
    dx = (end_x - start_x) / (N_x - 2);
}

Parameters::Parameters(double _density, double _g, double _Re, double _gamma, int _N_X, double _start_x, double _end_x, double _start_t, double _end_t, double _CFL, int _write_interval, int _max_iter, int _solver, int _initial, int _left_boundary, int _right_boundary) :
    density(_density),
    g(_g),
    Re(_Re),
    gamma(_gamma),
    N_x(_N_X),
    start_x(_start_x),
    end_x(_end_x),
    start_t(_start_t),
    end_t(_end_t),
    CFL(_CFL),
    write_interval(_write_interval),
    max_iter(_max_iter),
    solver(_solver),
    initial(_initial),
    left_boundary(_left_boundary),
    right_boundary(_right_boundary)
{
    c = 0.;
    calc_dx();
    calc_dt();
}

Parameters::Parameters() :
    density(0.),
    g(9.81),
    Re(0.),
    gamma(1.4),
    N_x(100),
    start_x(0.),
    end_x(1.),
    start_t(0.),
    end_t(1.),
    CFL(0.5),
    write_interval(5),
    max_iter(500),
    solver(1),
    initial(1),
    left_boundary(1),
    right_boundary(1)
{
    c = 0.;
    calc_dx();
    calc_dt();
}

Parameters::Parameters(std::string file_name) {
    std::ifstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << file_name << std::endl;
        Parameters();
        return;
    }

    std::unordered_map<std::string, std::string> params;
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, ' ') && std::getline(iss, value, ';')) {
            params[key] = value;
        }
    }
    file.close();

    //                                                       
    density = std::stod(params["density"]);
    g = std::stod(params["g"]);
    Re = std::stod(params["Re"]);
    gamma = std::stod(params["gamma"]);
    N_x = std::stoi(params["N_x"]) + 2; 
    start_x = std::stod(params["start_x"]);
    end_x = std::stod(params["end_x"]);
    start_t = std::stod(params["start_t"]);
    end_t = std::stod(params["end_t"]);
    CFL = std::stod(params["CFL"]);
    write_interval = std::stoi(params["writeInterval"]);
    max_iter = std::stoi(params["max_iter"]);
    solver = std::stoi(params["solver"]);
    initial = std::stoi(params["initial"]);
    left_boundary = std::stoi(params["left_boundary"]);
    right_boundary = std::stoi(params["right_boundary"]);
    c = 0.;
    calc_dx();
    calc_dt();
}

Parameters::~Parameters()
{
}

void Parameters::print() const {
    std::cout << "density: " << density << std::endl;
    std::cout << "g: " << g << std::endl;
    std::cout << "Re: " << Re << std::endl;
    std::cout << "gamma: " << gamma << std::endl;
    std::cout << "N_x: " << N_x << std::endl;
    std::cout << "start_x: " << start_x << std::endl;
    std::cout << "end_x: " << end_x << std::endl;
    std::cout << "dx: " << dx << std::endl;
    std::cout << "start_t: " << start_t << std::endl;
    std::cout << "end_t: " << end_t << std::endl;
    std::cout << "dt: " << dt << std::endl;
    std::cout << "CFL: " << CFL << std::endl;
    std::cout << "write_interval: " << write_interval << std::endl;
    std::cout << "initial: " << initial << std::endl;
    std::cout << "left_boundary: " << left_boundary << std::endl;
    std::cout << "right_boundary: " << right_boundary << std::endl;
    std::cout << std::endl;
}