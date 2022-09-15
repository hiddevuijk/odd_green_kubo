#ifndef GUARD_Time_CORRELATION_H
#define GUARD_Time_CORRELATION_H

/*
  Sample three different kinds of correlation functions:

  Where A and B can be of any type (they must be the same type),
  that have * and + defined, and A * B and A + B
  must be of the same type as A and B.

  The three classes are derived from the base class "TimeCorrelation"

  The constructor takes a single unsigned int corresponding to the number 
  of time differences (number_of_time_steps_),
  that is, the number of elements in the correlation function.

  A sample can be taken by the method Sample(T A, T B).
  It is assumed that samples are taken on equally spaced time points.

  If the number of samples is less than number_of_time_steps_, the
  correlation function is empty.

  Example:
  Initialize time correlation of two objects of type T:
    TimeCorrelation1<T> c(number_of_time_steps);
  Take a sample:
    c.Sample(a,b); 
  where a and b are of template type T
   
*/

#include <vector>

template <class T>
class TimeCorrelation {
 public:
  TimeCorrelation(unsigned int number_of_time_steps,
                  T (*function)(T,T,T,T));

  void Sample(T A, T B);

  // at least one sample has been taken.
  // if it returns false, GetTimeCorrelationFunction returns
  // a vector with 0s.
  bool InitStageDone() const { return number_of_samples_ > 0;}

  std::vector<T> GetTimeCorrelationFunction() const;

 protected:  

  void SampleCABFromList();

  unsigned int number_of_time_steps_;
  int number_of_samples_;


  // A_list_ constains
  // [A(t+n+1), ... ,A(t+tau-1),A(t+tau),A(t),A(t+1),...,A(t+n)] 
  // the integer "list_last_index_" is the location of the last added A:
  // A(t+tau) = A_list_[list_last_index_]
  std::vector<T> A_list_;
  std::vector<T> B_list_;

  unsigned int list_last_index_;

  // vector with samlples of  the correlation function
  std::vector<T> c_AB_;
 
  T (*function)(T A0, T At, T B0, T Bt); 
};


//////////////////////
// Member functions //
//////////////////////
template <class T>
TimeCorrelation<T>::TimeCorrelation(unsigned int number_of_time_steps,
                                    T (*function)(T,T,T,T) )
  : number_of_time_steps_(number_of_time_steps),
    number_of_samples_(-1 * number_of_time_steps),
    A_list_(number_of_time_steps),
    B_list_(number_of_time_steps),
    list_last_index_(-1),
    c_AB_(number_of_time_steps_),
    function(function)
{} 

template <class T>
void TimeCorrelation<T>::Sample(T A, T B)
{
  // if A_list_ and B_list_ not filled yet
  if (number_of_samples_ < 0) {
    ++list_last_index_;
    A_list_[list_last_index_] = A;
    B_list_[list_last_index_] = B;
    ++number_of_samples_;
  } else {
    list_last_index_ = (list_last_index_ + 1) % number_of_time_steps_;
    A_list_[list_last_index_] = A; 
    B_list_[list_last_index_] = B; 
    SampleCABFromList();
  } 
}

template <class T>
std::vector<T> TimeCorrelation<T>::GetTimeCorrelationFunction() const
{

  std::vector<T> c_AB_temp = c_AB_;

  if (number_of_samples_ > 0) {
    for (unsigned int i = 0; i < number_of_time_steps_; ++i) {
      c_AB_temp[i] /= number_of_samples_;
    }
  }

  return c_AB_temp;
}

template <class T>
void TimeCorrelation<T>::SampleCABFromList()
{
  ++number_of_samples_;
  int delta_ti;
  for (unsigned int i = 0;
       i < number_of_time_steps_; ++i)
  {
    delta_ti = list_last_index_ - i;
    if ( delta_ti < 0) {
      delta_ti += number_of_time_steps_;
    }
    c_AB_[i] += function( A_list_[list_last_index_],
                          A_list_[delta_ti],
                          B_list_[list_last_index_],
                          B_list_[delta_ti]);
  }
}


#endif
