
#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <boost/random.hpp>


class System {
 public:
	System(unsigned long int seed, double dt);

  void MakeTimeStep();
  void Integrate(int Ndt);

  double GetPosition() const { return x_;}
  double GetVelocity() const { return v_;}


  void ResetTime() { time_ = 0; }
  double GetTime() const { return time_;}

 private:
  const boost::normal_distribution<double> normal_distribution_;

  boost::mt19937 random_number_generator_;
  boost::variate_generator<boost::mt19937&,
        boost::normal_distribution<double> > random_normal_distribution_;

  // time step size
  double dt_;


	// particle velocity and position
  double v_, x_;
  double dx_t_, dx_mt_;


  // current time	
  double time_;

  double sqrt_2dt_;
};

System::System(unsigned long int seed, double dt)
  : normal_distribution_(0.0,1.0),
	  random_number_generator_(seed),
    random_normal_distribution_(random_number_generator_,
                                normal_distribution_),
	  dt_(dt), v_(0.0), x_(0.0), dx_t_(0.0), dx_mt_(0.0),
    time_(0.0), sqrt_2dt_(sqrt(2*dt))
{
  dx_mt_ = sqrt_2dt_ * random_normal_distribution_();
  dx_t_ = sqrt_2dt_ * random_normal_distribution_();
}


void System::Integrate(int Ndt)
{
  for (int i = 0; i < Ndt; ++i) MakeTimeStep();
}

void System::MakeTimeStep()
{
  x_ += dx_mt_; // x_ = x(t)

  dx_mt_ = dx_t_;
  dx_t_ = sqrt_2dt_ * random_normal_distribution_();

  v_ = (dx_t_ + dx_mt_) / (2 * dt_); // v_ = v(t)

  time_ += dt_;
}

#endif
