#ifndef GUARD_LOG_TIME_CORRELATION_H
#define GUARD_LOG_TIME_CORRELATION_H

#include "time_correlation.h"
#include <cmath>

template <class T>
class LogTimeCorrelation
{
 public:
  LogTimeCorrelation(unsigned int n_time_steps,
                     unsigned int number_of_decades,
                     T (*function)(T,T,T,T));

  void Sample(T A, T B);

  std::vector<T> GetTimeCorrelationFunction() const;
  std::vector<double> GetTimeList() const;
 private:
    unsigned int number_of_time_steps_;
    unsigned int number_of_decades_;

    std::vector<TimeCorrelation<T> > c_AB_list_;

    unsigned int number_of_samples_;
};

template <class T>
LogTimeCorrelation<T>::LogTimeCorrelation(unsigned int n_time_steps,
                                          unsigned int number_of_decades,
                                          T (*function)(T,T,T,T))
  : number_of_time_steps_(n_time_steps),
    number_of_decades_(number_of_decades),
    c_AB_list_(number_of_decades,
    TimeCorrelation<T>(n_time_steps,function)),
    number_of_samples_(0)
{}


template <class T>
void LogTimeCorrelation<T>::Sample(T A, T B)
{

  unsigned int D = 1;
  for (unsigned int dec = 0; dec < number_of_decades_; ++dec) {
    if (number_of_samples_ % D == 0) {
      c_AB_list_[dec].Sample(A, B);
    }
    D *= number_of_time_steps_;
  }
  number_of_samples_ += 1;
}

template <class T>
std::vector<T> LogTimeCorrelation<T>::GetTimeCorrelationFunction() const
{
   std::vector<T> c_AB_;
   for (unsigned int dec = 0; dec < number_of_decades_; ++dec) {
     std::vector<T> c_AB_decade =
           c_AB_list_[dec].GetTimeCorrelationFunction();
     for(unsigned int i = 0; i < number_of_time_steps_; ++i) {
       if ( (i == 0 and dec == 0) or i != 0) {
         c_AB_.push_back(c_AB_decade[i]);
       } 
     }
   }
  return c_AB_;
}

template <class T>
std::vector<double> LogTimeCorrelation<T>::GetTimeList() const
{

   std::vector<double> time_list;
   time_list.push_back(0.0);

   for (unsigned int dec = 0; dec < number_of_decades_; ++dec) {
     for(unsigned int i = 1; i < number_of_time_steps_; ++i) {
       time_list.push_back( i * pow(number_of_time_steps_,dec) );
     } 
   }

   return time_list;
}

#endif
