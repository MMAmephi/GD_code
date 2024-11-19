#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

class Parameters
{
public:
    double density; //���������
    double g;   //��������� ���������� �������
    double Re;  //����� ����������
    double gamma;   //���������� ��������
    double c; //�������� �����

    int N_x;
    double start_x;
    double end_x;
    double dx;

    double start_t;
    double end_t;
    double dt;

    double CFL;

    int write_interval;
    int max_iter;

    int solver;
    int initial;    //��� ���������� �������
    int left_boundary;  //��� ���������� ������� �����
    int right_boundary;  //��� ���������� ������� ������

    Parameters();

    Parameters(double _density, double _g, double _Re, double _gamma, int _N_X, double _start_x, double _end_x, double _start_t, double _end_t, double _CFL, int _write_interval, int _max_iter, int _solver, int _initial, int _left_boundary, int _right_boundary);

    Parameters(std::string file_name);

    void calc_dt();
    void calc_dx();

    void print() const;

    ~Parameters();
};