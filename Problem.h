#pragma once
#include "Parameters.h"
#include <filesystem>

namespace fs = std::filesystem;
typedef std::vector<double> vec;

class Problem
{
public:
	Parameters params;

	vec x;
	vec xc;

	vec x_new;

	vec rho;
	vec v_x;
	vec p;
	vec pb;

	vec rho_new;
	vec v_x_new;
	vec p_new;
	vec pb_new;


	vec m;
	vec imp_x;
	vec e;
	vec ei;

	vec imp_x_new;
	vec e_new;
	vec ei_new;

	vec w;

	Problem();

	Problem(Parameters& _params);

	~Problem();

	void write_out(int _n, double _time);

	void init_grid();
	void init();
	void calc_xc();
	void get_dt();

	void Test1();
	void Test2();
	void Test3();
	void Test4();
	void Test5();

	void Solve();
	void CrossSolverNonDivergent();
	void CrossSolverDivergent();
	void ArtificialViscosity();
	void ArtificialViscosityConst();
	void ArtificialViscosityLinear();
	void ArtificialViscosityNeumann();
	void ArtificialViscosityQuadro();
	void ArtificialViscosityMixed();
	//void Boundary_Wall(int _n);
};