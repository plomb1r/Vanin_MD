#include <iostream>
#include <string>
#include <fstream>
#include "params.h"
#include "constants.h"
#include "global_var.h"
#include "start_cond.h"
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <vector>
using std::vector;
using std::string;
using std::setprecision;

void getStrings() {
    string empty;
    std::ofstream StreamOut;
    StreamOut.open("G:\\Vanin_MD_2.txt");
    if (StreamOut.is_open())
    {
         StreamOut << empty << "\n"
            << "Number of Particles = " << NUMBERPARTICLES << "\n"
            << "Number of Steps =  " << NSTEPS << "\n"
            << "Size of System on X =  " << LX << "\n"
            << "Boltzmann Constant  =  " << setprecision(10) << K_B << "\n"
            << "Volume  =  " << VOLUME;
    }
}

vector<double> CalculationVectors() { //вычисление векторов и модулей rx12, ry12, rz12, rx21, ry21, rz21, r12_abs, r21_abs
    double rx12 = 0.0;
    double ry12 = 0.0;
    double rz12 = 0.0;
    double rx21 = 0.0;
    double ry21 = 0.0;
    double rz21 = 0.0;
    double r12_abs = 0.0;
    double r21_abs = 0.0;
    vector<double> resultVector;
    rx12 = coordx[0] - coordx[1];
    ry12 = coordy[0] - coordy[1];
    rz12 = coordz[0] - coordz[1];
    rx21 = coordx[1] - coordx[0];
    ry21 = coordy[1] - coordy[0];
    rz21 = coordz[1] - coordz[0];
    r12_abs = sqrt(rx12 * rx12 + ry12 * ry12 + rz12 * rz12);  //ƒлина вектора r12
    r21_abs = sqrt(rx21 * rx21 + ry21 * ry21 + rz21 * rz21);  //ƒлина вектора r21
    resultVector = { rx12, ry12, rz12, rx21, ry21, rz21, r12_abs, r21_abs }; //резултирующие данные по всем векторам
    return resultVector;

}

void PrintVectorsData(vector<double> v) {
    std::ofstream StreamOut;
    StreamOut.open("G:\\Vanin_MD_3.txt");
    if (StreamOut.is_open())
    {
        StreamOut << "Step = 0" << "\n"
            << "r1 = (rx1; ry1; rz1) = " << "(" << coordx[0] << "; " << coordy[0] << "; " << coordz[0] << ")" << "\n"
            << "r2 = (rx2; ry2; rz2) = " << "(" << coordx[1] << "; " << coordy[1] << "; " << coordz[1] << ")" << "\n"
            << "r12 = (rx12; ry12; rz12) = " << "(" << v[0] << "; " << v[1] << "; " << v[2] << ")" << "\n"
            << "r21 = (rx21; ry21; rz21) = " << "(" << v[3] << "; " << v[4] << "; " << v[5] << ")" << "\n"
            << "r12_abs = "  << v[6] << "\n"
            << "r21_abs = "  << v[7] << "\n"
            << "(rx21; ry21; rz21) / r12_abs = " << "(" << v[3]/v[6] << "; " << v[4]/v[6] << "; " << v[4]/v[6] << ")" << "\n"
            << "(rx12; ry12; rz12) / r12_abs = " << "(" << v[0]/v[6] << "; " << v[1]/v[6] << "; " << v[2]/v[6] << ")";
    }

}
void MD() { //вызов функции Molecular Dynamics, в которой создаетс€ и декларируетс€ экземпл€р класса vector, вызываетс€ функци€ н.у. дл€ 2 частиц
    vector<double> v;
    start_cond_two_particles();
    v = CalculationVectors(); //вызов метода подсчета данных дл€ вектора
    PrintVectorsData(v);
    //getStrings();
}

int main()
{
    MD();
    return 0;
}




