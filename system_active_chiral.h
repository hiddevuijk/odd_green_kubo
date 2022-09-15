
#ifndef GUARD_SYSTEM_ACTIVE_CHIRAL_H
#define GUARD_SYSTEM_ACTIVE_CHIRAL_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <boost/random.hpp>


class ActiveChiral {
 public:
	ActiveChiral(unsigned long int seed, double dt, double vs,
              double D, double Dr, double omega);

  void MakeTimeStep(double dt);
  void Integrate(double delta_t);

  double GetXPosition() const { return x_;}
  double GetYPosition() const { return y_;}
  double GetXVelocity() const { return vx_;}
  double GetYVelocity() const { return vy_;}


  void ResetTime() { time_ = 0; }
  double GetTime() const { return time_;}

 private:
  const boost::normal_distribution<double> normal_distribution_;

  boost::mt19937 random_number_generator_;
  boost::variate_generator<boost::mt19937&,
        boost::normal_distribution<double> > random_normal_distribution_;

  // time step size
  double dt_;

  double vs_, D_, Dr_, omega_;

	// particle velocity and position
  double x_, y_, vx_, vy_, dx_, dy_;
  double theta_;

  // current time	
  double time_;
};

ActiveChiral::ActiveChiral(unsigned long int seed, double dt,
               double vs, double D, double Dr, double omega)
  : normal_distribution_(0.0,1.0),
	  random_number_generator_(seed),
    random_normal_distribution_(random_number_generator_,
                                normal_distribution_),
	  dt_(dt), vs_(vs), D_(D), Dr_(Dr), omega_(omega),
    x_(0.0), y_(0.0), vx_(0.0), vy_(0.0), dx_(0.0), dy_(0.0),
    theta_(0.0), time_(0.0)
{}


void ActiveChiral::Integrate(double delta_time)
{
  while (delta_time > dt_) {
    MakeTimeStep(dt_);
    delta_time -= dt_;
  }
  MakeTimeStep(delta_time);
}

void ActiveChiral::MakeTimeStep(double dt)
{
  double sqrt_2_D_dt = sqrt(2 * D_ * dt);

  vx_ = dx_;
  vy_ = dy_;

  dx_ = vs_ * std::cos(theta_) * dt;
  dx_ += sqrt_2_D_dt * random_normal_distribution_();

  dy_ = vs_ * std::sin(theta_) * dt;
  dy_ += sqrt_2_D_dt *  random_normal_distribution_();

  vx_ += dx_;
  vx_ /= 2 * dt;
  vy_ += dy_;
  vy_ /= 2 * dt;

  x_ += dx_;
  y_ += dy_;

  theta_ += omega_ * dt +
            sqrt(2 * Dr_ * dt) * random_normal_distribution_();

  time_ += dt;
}

#endif
