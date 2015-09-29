#pragma once

#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cassert>

#include "definitions.hpp"

namespace consus
{

/// Discretize a number line from [Axis_Min_;Axis_Max_+0.5*Axis_Step_]
///
/// values corresponds to the x-value
/// entry corresponds to the y-value
template <class T = double>
class DiscreteAxis
{
  public:
   typedef typename std::vector<T>::const_iterator const_iterator;
   typedef typename std::vector<T>::iterator iterator;

   public:
    /// Minimal value (not center of first bin)
    T Min_;
    /// bin width
    T Step_;
    /// inverse Step_ width
    T Inv_Step_;
    /// number of
    int Num_Bins_;
    /// Max_imal value (not center of last bin)
    T Max_;
    /// internal array, stores a number for each bin
    std::vector<T> Array_;

  public:
   /// construct an Array over a discretized range [Min, Max + e[
   ///
   /// \param Min - first bin starts with this value
   /// \param Max - the upper bound will be the next possible value
   ///               greater or euqal then Max
   /// \param Step - bin width of every bin
   DiscreteAxis(const T Center_Min, const T Center_Max, const T Step,
                const T Initial_Value = 0.0)
       : Min_(Center_Min - 0.5 * Step),
         Step_(Step),
         Inv_Step_(1.0 / Step),
         Num_Bins_(calc_num_bins(Min_, Step, Center_Max)),
         Max_(Min_ + (Num_Bins_) * Step_),
         Array_(Num_Bins_, Initial_Value) {
      // test if Max value is in accessible range
     assert(Min_ < Max_);
     assert(Num_Bins_ > 0);
     assert(this->get_bin(this->get_max() - Step_/100.0) < Num_Bins_);
     
   };
    /// empty ctor
    DiscreteAxis(){};

    /// create default copy, move and assign constructors
    DiscreteAxis(const DiscreteAxis&) = default;
    DiscreteAxis(DiscreteAxis&&) = default;
    DiscreteAxis& operator=(const DiscreteAxis&) = default;
    DiscreteAxis& operator=(DiscreteAxis&&) = default;

    /// return the array location(array interval) for a given T value
    int get_bin(const T value) const noexcept {
      return static_cast<int>((value - Min_) * Inv_Step_);
    };

    /// return iterator pointing to the bin of a given value
    const_iterator iterator_at(const T value) const noexcept {
      return Array_.cbegin() + static_cast<int>((value - Min_) * Inv_Step_);
    }

    /* get functions in order to check ranges*/
    /// returns the lower bound of the first bin
    T get_min() const noexcept {
      return Min_;
    };

    /// returns the upper bound of the last bin
    T get_max() const noexcept {
      return Max_;
    };

    T get_center_min() const noexcept{
      return Min_ + 0.5 * Step_;
    }

    T get_center_max() const noexcept{
      return Max_ - 0.5 * Step_;
    }

    /// returns the bin width of every bin
    T get_step() const noexcept {
      return Step_;
    };

    /// returns number of bin
    int get_num_bins() const noexcept {
      return Num_Bins_;
    };

    /// return the center value of bin bin_
    T get_value(const int bin) const noexcept {
      return Min_ + Step_ * (bin + 0.5);
    };

    /// returns the entry for value value_
    T get_entry(const T value) const noexcept {
      return Array_[get_bin(value)];
    };

    /// add 1 to the corresponding bin
    void add_one(const T energy) noexcept {
      this->Array_[get_bin(energy)] += 1.0;
    };
    
    /// reset histogram entries to zero;
    /// replace with std algorithm
    void reset(const T value = 0.0) noexcept {
      std::fill(Array_.begin(), Array_.end(), value);
    };

    /* access Array_ elements directly with the [] operator*/
    /// get Referennce to object in order to access const element
    const T& operator [] (const int entry) const {
      return Array_[entry];
    }
    /// get Reference to object in order to read and write elements
    T& operator [] (const int entry) {
      return Array_[entry];
    }

    /// returns pointer to internal data
    T *data() noexcept {
      return Array_.data();
    };
    /// returns const pointer to internal data
    const T* data() const noexcept {
      return Array_.data();
    };

    /// returns the number of elements in the object
    int size() const noexcept { return Array_.size(); }
    /// returns the Iterator to the beginning of DiscreteAxis
    iterator begin() { return Array_.begin(); }
    /// returns the Iteratot to the end of DiscreteAxis
    iterator end() { return Array_.end(); }
    /// const_iterator to first element of internal DiscreteAxis
    const_iterator begin() const { return Array_.begin(); }
    /// const_iterator to one after last element of internal DiscreteAxis
    const_iterator end() const { return Array_.end(); }
    /// const_iterator to first element of internal DiscreteAxis
    const_iterator cbegin() const { return Array_.cbegin(); }
    /// const_iterator to one after last element of internal DiscreteAxis
    const_iterator cend() const { return Array_.cend(); }

    /// calculate number of bins
    int calc_num_bins(const T Min, const T Step,
                      const T Max) const;

    /*-----------------------------------------------------------------------------
     *  comparision operator
     *-----------------------------------------------------------------------------*/
    /// test for equality
    bool operator==(const DiscreteAxis& rhs) const noexcept;

    /// test for inequality
    bool operator!=(const DiscreteAxis& rhs) const noexcept;

    /// returns a formatted string
    std::string formatted_string() const;

    /// write the object into a file such that gnuplot can read it
    void write_formatted(const std::string& file) const;
    /// load a file which is formatted according to write_formatted
    void load_formatted(const std::string& file);
};

template <class T>
inline bool DiscreteAxis<T>::operator==(const DiscreteAxis<T>& rhs) const noexcept{
  if (Min_ != rhs.Min_) {
    return false;
  }
  if (Step_ != rhs.Step_) {
    return false;
  }
  if (Num_Bins_ != rhs.Num_Bins_) {
    return false;
  }
  if (Max_ != rhs.Max_) {
    return false;
  }
  if (Array_ != rhs.Array_) {
    return false;
  }
  return true;
}

template <class T>
inline bool DiscreteAxis<T>::operator!=(const DiscreteAxis& rhs) const noexcept{
  return !(*this == rhs);
}

template <class T>
inline int DiscreteAxis<T>::calc_num_bins(const T Min, const T Step,
                                          const T Max) const {
  assert(Min < Max);
  T num = (Max-Min)/Step;
  return num + 1;
}

template <class T>
inline std::string DiscreteAxis<T>::formatted_string() const {
  std::stringstream ss;
  ss.precision(18);
  ss.setf(std::ios::scientific);
  ss << "#AxisMin= "  << get_min() << "\n";
  ss << "#AxisMax= "  << get_max() << "\n";
  ss << "#AxisStep= " << Step_ << "\n";
  ss << "#NumBins= " << Num_Bins_ << "\n";
  for (int i = 0; i < Num_Bins_; i++) {
    ss << get_value(i) << " " << Array_[i] << "\n";
  }
  return ss.str();
}

template <class T>
inline void DiscreteAxis<T>::write_formatted(const std::string& filename) const {
  std::fstream out(filename, std::ios::out | std::ios::trunc);
  out << this->formatted_string();
}

template <class T>
void  DiscreteAxis<T>::load_formatted(const std::string& filename)
{
  std::fstream file(filename, std::ios::in);
  std::string name;
  file >> name;
  if (name != "#AxisMin="){
    std::cerr << "Error while load_formatted of DiscreteAxis " << filename
              << "\n";
    std::cerr << "#AxisMin= is missing\n";
    std::exit(1);
  }
  file >> Min_;
  file >> name;
  if (name != "#AxisMax="){
    std::cerr << "Error while load_formatted of DiscreteAxis " << filename
              << "\n";
    std::cerr << "#AxisMax= is missing\n";
    std::exit(1);
  }
  file >> Max_;
  file >> name;
  if (name != "#AxisStep="){
    std::cerr << "Error while load_formatted of DiscreteAxis " << filename
              << "\n";
    std::cerr << "#AxisStep= is missing\n";
    std::exit(1);
  }
  file >> Step_;
  Inv_Step_ = 1.0/Step_;
  file >> name;
  if (name != "#NumBins="){
    std::cerr << "Error while load_formatted of DiscreteAxis " << filename
              << "\n";
    std::cerr << "#NumBins= is missing\n";
    std::exit(1);
  }
  file >> Num_Bins_;
  Array_.resize(Num_Bins_);
  size_t lines = 0;
  T value_tmp_old = 0.0;
  std::string line;
  while(std::getline(file, line)){
    T value_tmp_new = 0.0, entry_tmp = 0.0;
    // test if only valid chars are in the line
    if (line[0] != '#' and not line.empty()) {
      std::stringstream(line) >> value_tmp_new >> entry_tmp;
      if(lines > Array_.size()-1){
        std::cerr << "Error while load_formatted of DiscreteAxis " << filename
                  << "\n";
        std::cerr << "Number of lines does not fit to NumBins\n";
        std::cerr << Num_Bins_ << " " << lines << " " << Array_.size()-1 << "\n";
        std::exit(1);
      }
      Array_[lines] = entry_tmp;
      if (lines >= 1){
        T Step_tmp = value_tmp_new - value_tmp_old;
        if (std::fabs(Step_ - Step_tmp) > 1e-7) {
          std::cerr << "Error while load_formatted of DiscreteAxis " << filename
                    << "\n";
          std::cerr << "Step_ size is not consistent in line " << lines << "\n";
          std::cerr << value_tmp_new << " " << value_tmp_old << "\n";
          std::cerr << Step_ << " vs. " << Step_tmp << "\n";
        }
      }
      value_tmp_old = value_tmp_new;
      ++lines;
    }
  }
}

template <class T>
inline void serialize(const DiscreteAxis<T>& t, std::ostream& out){
    serialize(t.Min_, out);
    serialize(t.Step_, out);
    serialize(t.Inv_Step_, out);
    serialize(t.Num_Bins_, out);
    serialize(t.Max_, out);
    serialize(t.Array_, out);
}

template <class T>
inline void deserialize(DiscreteAxis<T>& t, std::istream& in){
    deserialize(t.Min_, in);
    deserialize(t.Step_, in);
    deserialize(t.Inv_Step_, in);
    deserialize(t.Num_Bins_, in);
    deserialize(t.Max_, in);
    deserialize(t.Array_, in);
}

} /* end of namespace betamc */ 
