
#include "time_correlation.h"
#include "log_time_correlation.h"
#include "system_active_chiral.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>


using namespace std;

template <class T>
void save_vector(const vector<T>& vec, string name)
{
  ofstream vec_out; 
  vec_out.open(name);

  for (unsigned int i = 0; i < vec.size(); ++i) {
    vec_out << vec[i] << "\n";
  }

  vec_out.close();
}

double frr(double A0, double At, double B0, double Bt)
{
  return (At - A0) * (Bt - B0);
}

double fvr(double A0, double At, double B0, double Bt)
{
  return A0 * (B0 - Bt);
}

double fvv(double A0, double At, double B0, double Bt)
{
  return A0 * Bt;
}

int main()
{
  double dt = 1e-2;
  double t_sample = dt;
  double t_total = 1e6;

  unsigned int n_time_steps = 5;
  unsigned int n_decades = 5;

  long unsigned int seed = 213451789;

  double vs = 1.0;
  double D = 0.5;
  double Dr = 1.0;
  double omega = 1.0;

  ActiveChiral active_chiral(seed, dt, vs, D, Dr, omega);

  LogTimeCorrelation<double> cxx(n_time_steps, n_decades, frr); 
  LogTimeCorrelation<double> cyy(n_time_steps, n_decades, frr); 
  LogTimeCorrelation<double> cxy(n_time_steps, n_decades, frr); 
  LogTimeCorrelation<double> cyx(n_time_steps, n_decades, frr); 
  LogTimeCorrelation<double> cvxx(n_time_steps, n_decades, fvr); 
  LogTimeCorrelation<double> cvyy(n_time_steps, n_decades, fvr); 
  LogTimeCorrelation<double> cvxy(n_time_steps, n_decades, fvr); 
  LogTimeCorrelation<double> cvyx(n_time_steps, n_decades, fvr); 
  LogTimeCorrelation<double> cvxvx(n_time_steps, n_decades, fvv); 
  LogTimeCorrelation<double> cvyvy(n_time_steps, n_decades, fvv); 
  LogTimeCorrelation<double> cvxvy(n_time_steps, n_decades, fvv); 
  LogTimeCorrelation<double> cvyvx(n_time_steps, n_decades, fvv); 
 

  active_chiral.Integrate(t_sample);
  while (active_chiral.GetTime() < t_total) {
    active_chiral.Integrate(t_sample);

    cxx.Sample(active_chiral.GetXPosition(),
               active_chiral.GetXPosition()); 

    cyy.Sample(active_chiral.GetYPosition(),
               active_chiral.GetYPosition()); 

    cxy.Sample(active_chiral.GetXPosition(),
               active_chiral.GetYPosition()); 

    cyx.Sample(active_chiral.GetYPosition(),
               active_chiral.GetXPosition()); 

    cvxx.Sample(active_chiral.GetXVelocity(),
                active_chiral.GetXPosition()); 

    cvyy.Sample(active_chiral.GetYVelocity(),
                active_chiral.GetYPosition()); 

    cvxy.Sample(active_chiral.GetXVelocity(),
                active_chiral.GetYPosition()); 

    cvyx.Sample(active_chiral.GetYVelocity(),
                active_chiral.GetXPosition()); 

    cvxvx.Sample(active_chiral.GetXVelocity(),
                 active_chiral.GetXVelocity()); 

    cvyvy.Sample(active_chiral.GetYVelocity(),
                 active_chiral.GetYVelocity()); 

    cvxvy.Sample(active_chiral.GetXVelocity(),
                 active_chiral.GetYVelocity()); 

    cvyvx.Sample(active_chiral.GetYVelocity(),
                 active_chiral.GetXVelocity()); 

  
    cout << t_total << "\t" << active_chiral.GetTime() << endl;
  } 

  save_vector(cxx.GetTimeList(), "t.dat");
  save_vector(cxx.GetTimeCorrelationFunction(), "cxx.dat");
  save_vector(cyy.GetTimeCorrelationFunction(), "cyy.dat");
  save_vector(cxy.GetTimeCorrelationFunction(), "cxy.dat");
  save_vector(cyx.GetTimeCorrelationFunction(), "cyx.dat");
  save_vector(cvxx.GetTimeCorrelationFunction(), "cvxx.dat");
  save_vector(cvyy.GetTimeCorrelationFunction(), "cvyy.dat");
  save_vector(cvxy.GetTimeCorrelationFunction(), "cvxy.dat");
  save_vector(cvyx.GetTimeCorrelationFunction(), "cvyx.dat");
  save_vector(cvxvx.GetTimeCorrelationFunction(), "cvxvx.dat");
  save_vector(cvyvy.GetTimeCorrelationFunction(), "cvyvy.dat");
  save_vector(cvxvy.GetTimeCorrelationFunction(), "cvxvy.dat");
  save_vector(cvyvx.GetTimeCorrelationFunction(), "cvyvx.dat");


  return 0;
}
