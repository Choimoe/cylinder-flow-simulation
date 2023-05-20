#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdio.h>

double d1 = 2;
double d2 = 40;

const int M = 100;
const int N = 200;

const double PI = 3.1415926535897932384626;

MatrixXd X(M + 1, N + 1);
MatrixXd Y(M + 1, N + 1);

int main() {
  double theta;
  double t = 1.032688205;
  for (int i = 0; i <= M; i++) {
    theta = 2 * PI * float(i) / float(M);
    for (int j = 0; j <= N; j++) {
      double r = d1 / 2.0 + 0.001*(1.0 - pow(t, j)) / (1.0 - t);
      X(i, j) = r * cos(theta);
      Y(i, j) = r * sin(theta);
    }
  }

  ofstream Xout("X_new.txt");
  Xout << X;
  Xout.close();
  ofstream Yout("Y_new.txt");
  Yout << Y;
  Yout.close();
}
