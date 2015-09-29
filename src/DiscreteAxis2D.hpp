#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

#include "Range.hpp"
#include "definitions.hpp"

namespace consus
{
  
/// class for storing a 2 dimensional histogram
///
/// data is stored internally in a one dimensional container,
/// the second dimension is stored continuously
template <class T = double>
class DiscreteAxis2D {
public:
  typedef typename std::vector<T>::const_iterator const_iterator;
  typedef typename std::vector<T>::iterator iterator;

public:
  T Min_First_, Max_First_, Step_First_;
  T Center_Min_First_;
  T Inv_Step_First_;
  T Min_Second_, Max_Second_, Step_Second_;
  T Center_Min_Second_;
  T Inv_Step_Second_;
  int Num_Bins_First_, Num_Bins_Second_;
  std::vector<T> Array_;

public:
  /// constructor
  DiscreteAxis2D();
  DiscreteAxis2D(const T Center_Min_First, const T Center_Max_First,
                 const T Step_First, const T Center_Min_Second,
                 const T Center_Max_Second, const T Step_Second,
                 const T Initial_Value = 0);
  DiscreteAxis2D(const Range<T> &rangeE1, const Range<T> &rangeE2,
                 const T Initial_Value = 0)
      : DiscreteAxis2D(rangeE1.start, rangeE1.step, rangeE1.end, rangeE2.start,
                       rangeE2.step, rangeE2.end, Initial_Value){};

  /// return the array index for a set of values
  int get_bin(const T value_first, const T value_second) const;

  int get_bin_first(const T value) const {
    assert(value >= Min_First_);
    assert(value <= Max_First_ + Step_First_);
    return static_cast<int>((value - Min_First_) * Inv_Step_First_);
  }

  int get_bin_second(const T value) const {
    assert(value >= Min_Second_);
    assert(value <= Max_Second_ + Step_Second_);
    return static_cast<int>((value - Min_Second_) * Inv_Step_Second_);
  }

  /// see above
  int get_bin(const std::pair<T, T> &value) const;

  int get_index(const int index_first, const int index_second) const {
    return index_second + index_first * Num_Bins_Second_;
  }

  T get_entry(const int index_first, const int index_second) const {
    return Array_[get_index(index_first, index_second)];
  }

  T get_entry(const int index) const{
    return Array_[index];
  }

  /// increment operator on histogram element corresponding to value_first and
  /// value_second
  void add_one(const std::pair<T, T> &value) {
    Array_[get_bin(value)]++;
  }

  void add_one(const T value_first, const T value_second){
    Array_[get_bin(value_first, value_second)]++;
  }

  /// return minimal range of first dimension
  T get_min_first() const { return Min_First_; }
  /// return center of first bin of first dimension
  T get_center_min_first() const { return Center_Min_First_; }
  /// return minimal range of second dimension
  T get_min_second() const { return Min_Second_; }
  /// return center of first bin of second dimension
  T get_center_min_second() const { return Center_Min_Second_; }
  /// return maximal range of first dimension
  T get_max_first() const { return Max_First_; }
  /// return center of last bin of first dimension
  T get_center_max_first() const { return Max_First_ - 0.5 * Step_First_; }
  /// return maximal range of second dimension
  T get_max_second() const { return Max_Second_; }
  /// return center of last bin of second dimension
  T get_center_max_second() const { return Max_Second_ - 0.5 * Step_Second_; }
  /// get step range of first dimension
  T get_step_first() const { return Step_First_; }
  /// get step range of second dimension
  T get_step_second() const { return Step_Second_; }
  /// get number of bins in first dimension
  int get_num_bins_first() const { return Num_Bins_First_; }
  /// get number of bins in second dimension
  int get_num_bins_second() const { return Num_Bins_Second_; }
  /// get corresponding values of bin with index bin_
  std::pair<T, T> get_value(const int bin) const {
    int pos_first = bin / Num_Bins_First_;
    int pos_second = bin % Num_Bins_First_;
    return std::make_pair(Center_Min_First_ + Step_First_ * pos_first,
                          Center_Min_Second_ + Step_Second_ * pos_second);
  }
  /// get corresponding value from first dimension of bin  with index bin_
  T get_value_first(const int bin_first) const {
    return Center_Min_First_ + Step_First_ * bin_first;
  }
  /// get corresponding value from second dimension of bin  with index bin_
  T get_value_second(const int bin_second) const {
    return Center_Min_Second_ + Step_Second_ * bin_second;
  }

  Range<T> get_range_first() const {
    return Range<T>{Center_Min_First_, Step_First_, Max_First_-0.5*Step_First_};
  }

  Range<T> get_range_second() const {
    return Range<T>{Min_Second_, Step_Second_, Max_Second_};
  }

  bool is_in_boundaries(const T value_first, const T value_second){
    return (value_first >= Min_First_ and value_first <
            Max_First_ and value_second >= Min_Second_ and value_second <
            Max_Second_);
  }

  /// set every entry to 0
  void reset() { std::fill(Array_.begin(), Array_.end(), T(0)); }

  /*-----------------------------------------------------------------------------
   *  add container like interface
   *-----------------------------------------------------------------------------*/
  /// get Referennce to object in order to access const element
  const T &operator[](const int entry) const { return this->Array_[entry]; }
  /// get Reference to object in order to read and write elements
  T &operator[](const int entry) { return this->Array_[entry]; }
  /// returns pointer to internal data
  T *data() noexcept { return Array_.data(); };
  /// returns const pointer to internal data
  const T *data() const noexcept { return Array_.data(); };

  /// returns the number of elements in the object
  int size() const { return Array_.size(); }
  /// returns the Iterator to the beginning of DiscreteAxis2D
  iterator begin() { return Array_.begin(); }
  /// returns the Iteratot to the end of DiscreteAxis2D
  iterator end() { return Array_.end(); }
  /// const_iterator to first element of internal DiscreteAxis2D
  const_iterator begin() const { return Array_.begin(); }
  /// const_iterator to one after last element of internal DiscreteAxis2D
  const_iterator end() const { return Array_.end(); }

  /*-----------------------------------------------------------------------------
   *  friend declaration
   *-----------------------------------------------------------------------------*/
  /// overload ostream operator
  friend std::ostream &operator<<(std::ostream &out,
                                  const DiscreteAxis2D &hist) {
    out << "#CenterMinFirst " << hist.get_center_min_first() << "\n"
        << "#CenterMaxFirst " << hist.get_center_max_first() << "\n"
        << "#StepFirst " << hist.get_step_first() << "\n"
        << "#CenterMinSecond " << hist.get_center_min_second() << "\n"
        << "#CenterMaxSecond " << hist.get_center_max_second() << "\n"
        << "#StepSecond " << hist.get_step_first() << "\n";
    for (int i = 0; i < hist.get_num_bins_first(); i++) {
      for (int j = 0; j < hist.get_num_bins_second(); ++j) {
        out << hist.get_value_first(i) << " "
            << hist.get_value_second(j) << " "
            << hist[hist.get_index(i,j)] << "\n";
      }
      out << "\n";
    }
    return out;
  }

  void write_formatted(const std::string& filename){
    std::ofstream file(filename, std::ios::out | std::ios::trunc);
    file << *this;
  }
};

template <class T>
DiscreteAxis2D<T>::DiscreteAxis2D(
    const T Center_Min_First, const T Center_Max_First, const T Step_First,
    const T Center_Min_Second, const T Center_Max_Second,
    const T Step_Second, const T Initial_Value)
    : Min_First_(Center_Min_First - 0.5 + Step_First),
      Max_First_(Center_Max_First + 0.499 * Step_First),
      Step_First_(Step_First),
      Center_Min_First_(Center_Max_First),
      Inv_Step_First_(1.0 / Step_First_),
      Min_Second_(Center_Min_Second - 0.5 + Step_Second_),
      Max_Second_(Center_Max_Second + 0.499 * Step_Second_),
      Step_Second_(Step_Second),
      Center_Min_Second_(Min_Second_ + 0.5 * Step_Second_),
      Inv_Step_Second_(1.0 / Step_Second),
      Num_Bins_First_(
          static_cast<int>((Max_First_ - Min_First_) / Step_First_) + 1),
      Num_Bins_Second_(
          static_cast<int>((Max_Second_ - Min_Second_) / Step_Second_) + 1),
      Array_(Num_Bins_First_ * Num_Bins_Second_, Initial_Value) {
  assert(Min_First_ < Max_First_);
  assert(Min_Second_ < Max_Second_);

  Max_First_ = (Min_First_ + Num_Bins_First_ * Step_First_);
  Max_Second_ = (Min_Second_ + Num_Bins_Second_ * Step_Second_);
}

template <class T>
DiscreteAxis2D<T>::DiscreteAxis2D(){}

template <class T>
inline int DiscreteAxis2D<T>::get_bin(const T value_first,
                                      const T value_second) const {
  assert(value_first >= Min_First_);
  assert(value_first <= Max_First_ + Step_First_);
  assert(value_second >= Min_Second_);
  assert(value_second <= Max_Second_ + Step_Second_);
  int pos_first =
      static_cast<int>((value_first - Min_First_) * Inv_Step_First_);
  int pos_second =
      static_cast<int>((value_second - Min_Second_) * Inv_Step_Second_);
  return pos_second + pos_first * Num_Bins_Second_;
}

template <class T>
inline int DiscreteAxis2D<T>::get_bin(const std::pair<T, T> &value)
    const {
  assert(value.first >= Min_First_);
  assert(value.first <= Max_First_ + Step_First_);
  assert(value.second >= Min_Second_);
  assert(value.second <= Max_Second_ + Step_Second_);
  int pos_first =
      static_cast<int>((value.first - Min_First_) * Inv_Step_First_);
  int pos_second =
      static_cast<int>((value.second - Min_Second_) * Inv_Step_Second_);
  return pos_second + pos_first * Num_Bins_Second_;
}

template <class T>
inline void serialize(const DiscreteAxis2D<T>& t, std::ostream& out) {
  dlib::serialize(t.Min_First_, out);
  dlib::serialize(t.Max_First_, out);
  dlib::serialize(t.Step_First_, out);
  dlib::serialize(t.Center_Min_First_, out);
  dlib::serialize(t.Inv_Step_First_, out);
  dlib::serialize(t.Min_Second_, out);
  dlib::serialize(t.Max_Second_, out);
  dlib::serialize(t.Step_Second_, out);
  dlib::serialize(t.Center_Min_Second_, out);
  dlib::serialize(t.Inv_Step_Second_, out);
  dlib::serialize(t.Num_Bins_First_, out);
  dlib::serialize(t.Num_Bins_Second_, out);
  dlib::serialize(t.Array_, out);
}

template <class T>
inline void deserialize(DiscreteAxis2D<T>& t, std::istream& in) {
  dlib::deserialize(t.Min_First_, in);
  dlib::deserialize(t.Max_First_, in);
  dlib::deserialize(t.Step_First_, in);
  dlib::deserialize(t.Center_Min_First_, in);
  dlib::deserialize(t.Inv_Step_First_, in);
  dlib::deserialize(t.Min_Second_, in);
  dlib::deserialize(t.Max_Second_, in);
  dlib::deserialize(t.Step_Second_, in);
  dlib::deserialize(t.Center_Min_Second_, in);
  dlib::deserialize(t.Inv_Step_Second_, in);
  dlib::deserialize(t.Num_Bins_First_, in);
  dlib::deserialize(t.Num_Bins_Second_, in);
  dlib::deserialize(t.Array_, in);
}

template <class T = double>
DiscreteAxis2D<T> make_DiscreteAxis2d(
    const T Min_first, const T Max_First, const T Step_First,
    const T Min_Second, const T Max_Second,
    const T Step_Second, const T Initial_Value = 0) {
  return DiscreteAxis2D<T>(Min_first, Max_First, Step_First, Min_Second,
                           Max_Second, Step_Second, Initial_Value);
}

template <class T = double>
DiscreteAxis2D<T> make_DiscreteAxis2d(const Range<T> &rangeE1,
                                        const Range<T> &rangeE2,
                                        const T Initial_Value = 0) {
  return DiscreteAxis2D<T>(rangeE1, rangeE2, Initial_Value);
}

} /* end of namespace consus */ 
