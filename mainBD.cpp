
#include "time_correlation.h"
#include "log_time_correlation.h"
#include "systemBD.h"

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

int main()
{
  double dt = 1e-1;
  double t_total = 1e5;


  unsigned int n_time_steps = 1000;

  long unsigned int seed = 113456789;


  System system(seed, dt);
  
  TimeCorrelation1<double> cvv(n_time_steps); 
  TimeCorrelation2<double> cvr(n_time_steps); 
  TimeCorrelation3<double> crr(n_time_steps); 

  vector<double> x,v;

  system.Integrate(100);
  while (system.GetTime() < t_total) {
    system.MakeTimeStep();
    x.push_back(system.GetPosition());
    v.push_back(system.GetVelocity());
    cvv.Sample(system.GetVelocity(), system.GetVelocity());
    cvr.Sample(system.GetVelocity(), system.GetPosition());
    crr.Sample(system.GetPosition(), system.GetPosition());
    cout << t_total << "\t" << system.GetTime() << endl;
  } 

  save_vector(x, "x.dat");
  save_vector(v, "v.dat");

  save_vector(cvv.GetTimeCorrelationFunction(), "cvv.dat");
  save_vector(cvr.GetTimeCorrelationFunction(), "cvr.dat");
  save_vector(crr.GetTimeCorrelationFunction(), "crr.dat");

  return 0;
}
