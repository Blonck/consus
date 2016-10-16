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

#include "Range.hpp"
#include "definitions.hpp"

namespace consus
{

/// Discretize a number line from [Axis_min_;Axis_max_+0.5*Axis_step_]
///
/// values corresponds to the x-value
/// entry corresponds to the y-value
class DiscreteAxis
{
  public:
   typedef typename std::vector<double>::const_iterator const_iterator;
   typedef typename std::vector<double>::iterator iterator;

   public:
    /// minimal value (not center of first bin)
    double min_;
    /// bin width
    double step_;
    /// inverse step_ width
    double invstep_;
    /// number of
    int numBins_;
    /// maximal value (not center of last bin)
    double max_;
    /// internal array, stores a number for each bin
    std::vector<double> array_;

  public:
   /// construct an Array over a discretized range [min, max + e[
   ///
   /// \param min - first bin starts with this value
   /// \param max - the upper bound will be the next possible value
   ///               greater or euqal then max
   /// \param step - bin width of every bin
   DiscreteAxis(const double min, const double max, const double step,
                const double initialValue = 0.0)
       : min_(min),
         step_(step),
         invstep_(1.0 / step),
         numBins_(calc_num_bins(min_, step, max)),
         max_(min_ + (numBins_) * step_),
         array_(numBins_, initialValue) {
      // test if max value is in accessible range
     assert(min_ < max_);
     assert(numBins_ > 0);
     assert(this->get_bin(this->get_max() - step_/100.0) < numBins_);
     
   };
   DiscreteAxis(const Range& range, const double initialValue = 0.0)
     : DiscreteAxis(range.start, range.end, range.step, initialValue){};

   /// empty ctor
   DiscreteAxis(){};

    /// create default copy, move and assign constructors
    DiscreteAxis(const DiscreteAxis&) = default;
    DiscreteAxis(DiscreteAxis&&) = default;
    DiscreteAxis& operator=(const DiscreteAxis&) = default;
    DiscreteAxis& operator=(DiscreteAxis&&) = default;

    /// return the array location(array interval) for a given double value
    int get_bin(const double value) const noexcept {
      return static_cast<int>((value - min_) * invstep_);
    };

    /// return iterator pointing to the bin of a given value
    const_iterator iterator_at(const double value) const noexcept {
      return array_.cbegin() + static_cast<int>((value - min_) * invstep_);
    }

    /* get functions in order to check ranges*/
    /// returns the lower bound of the first bin
    double get_min() const noexcept {
      return min_;
    };

    /// returns the upper bound of the last bin
    double get_max() const noexcept {
      return max_;
    };

    double get_center_min() const noexcept{
      return min_ + 0.5 * step_;
    }

    double get_center_max() const noexcept{
      return max_ - 0.5 * step_;
    }

    /// returns the bin width of every bin
    double get_step() const noexcept {
      return step_;
    };

    /// returns number of bin
    int get_num_bins() const noexcept {
      return numBins_;
    };

    /// return the center value of bin bin_
    double get_value(const int bin) const noexcept {
      return min_ + step_ * (bin + 0.5);
    };

    /// returns the entry for value value_
    double get_entry(const double value) const noexcept {
      return array_[get_bin(value)];
    };

    /// add 1 to the corresponding bin
    void add_one(const double energy) noexcept {
      this->array_[get_bin(energy)] += 1.0;
    };
    
    /// reset histogram entries to zero;
    /// replace with std algorithm
    void reset(const double value = 0.0) noexcept {
      std::fill(array_.begin(), array_.end(), value);
    };

    /* access array_ elements directly with the [] operator*/
    /// get Referennce to object in order to access const element
    const double& operator [] (const int entry) const {
      return array_[entry];
    }
    /// get Reference to object in order to read and write elements
    double& operator [] (const int entry) {
      return array_[entry];
    }

    /// returns pointer to internal data
    double *data() noexcept {
      return array_.data();
    };
    /// returns const pointer to internal data
    const double* data() const noexcept {
      return array_.data();
    };

    /// returns the number of elements in the object
    int size() const noexcept { return array_.size(); }
    /// returns the Iterator to the beginning of DiscreteAxis
    iterator begin() { return array_.begin(); }
    /// returns the Iteratot to the end of DiscreteAxis
    iterator end() { return array_.end(); }
    /// const_iterator to first element of internal DiscreteAxis
    const_iterator begin() const { return array_.begin(); }
    /// const_iterator to one after last element of internal DiscreteAxis
    const_iterator end() const { return array_.end(); }
    /// const_iterator to first element of internal DiscreteAxis
    const_iterator cbegin() const { return array_.cbegin(); }
    /// const_iterator to one after last element of internal DiscreteAxis
    const_iterator cend() const { return array_.cend(); }

    /// calculate number of bins
    int calc_num_bins(const double min, const double step,
                      const double max) const;

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

inline bool DiscreteAxis::operator==(const DiscreteAxis& rhs) const noexcept {
  if (min_ != rhs.min_) {
    return false;
  }
  if (step_ != rhs.step_) {
    return false;
  }
  if (numBins_ != rhs.numBins_) {
    return false;
  }
  if (max_ != rhs.max_) {
    return false;
  }
  if (array_ != rhs.array_) {
    return false;
  }
  return true;
}

inline bool DiscreteAxis::operator!=(const DiscreteAxis& rhs) const noexcept {
  return !(*this == rhs);
}

inline int DiscreteAxis::calc_num_bins(const double min, const double step,
                                       const double max) const {
  assert(min < max);
  double num = (max-min)/step;
  return num + 1;
}

inline std::string DiscreteAxis::formatted_string() const {
  std::stringstream ss;
  ss.precision(18);
  ss.setf(std::ios::scientific);
  ss << "#AxisMin= "  << get_min() << "\n";
  ss << "#AxisMax= "  << get_max() << "\n";
  ss << "#AxisStep= " << step_ << "\n";
  ss << "#NumBins= " << numBins_ << "\n";
  for (int i = 0; i < numBins_; i++) {
    ss << get_value(i) << " " << array_[i] << "\n";
  }
  return ss.str();
}

inline void DiscreteAxis::write_formatted(const std::string& filename) const {
  std::fstream out(filename, std::ios::out | std::ios::trunc);
  out << this->formatted_string();
}

void  DiscreteAxis::load_formatted(const std::string& filename)
{
  std::fstream file(filename, std::ios::in);
  std::string name;
  file >> name;
  if (name != "#Axismin="){
    std::cerr << "Error while load_formatted of DiscreteAxis " << filename
              << "\n";
    std::cerr << "#Axismin= is missing\n";
    std::exit(1);
  }
  file >> min_;
  file >> name;
  if (name != "#Axismax="){
    std::cerr << "Error while load_formatted of DiscreteAxis " << filename
              << "\n";
    std::cerr << "#Axismax= is missing\n";
    std::exit(1);
  }
  file >> max_;
  file >> name;
  if (name != "#Axisstep="){
    std::cerr << "Error while load_formatted of DiscreteAxis " << filename
              << "\n";
    std::cerr << "#Axisstep= is missing\n";
    std::exit(1);
  }
  file >> step_;
  invstep_ = 1.0/step_;
  file >> name;
  if (name != "#NumBins="){
    std::cerr << "Error while load_formatted of DiscreteAxis " << filename
              << "\n";
    std::cerr << "#NumBins= is missing\n";
    std::exit(1);
  }
  file >> numBins_;
  array_.resize(numBins_);
  size_t lines = 0;
  double value_tmp_old = 0.0;
  std::string line;
  while(std::getline(file, line)){
    double value_tmp_new = 0.0, entry_tmp = 0.0;
    // test if only valid chars are in the line
    if (line[0] != '#' and not line.empty()) {
      std::stringstream(line) >> value_tmp_new >> entry_tmp;
      if(lines > array_.size()-1){
        std::cerr << "Error while load_formatted of DiscreteAxis " << filename
                  << "\n";
        std::cerr << "Number of lines does not fit to NumBins\n";
        std::cerr << numBins_ << " " << lines << " " << array_.size()-1 << "\n";
        std::exit(1);
      }
      array_[lines] = entry_tmp;
      if (lines >= 1){
        double step_tmp = value_tmp_new - value_tmp_old;
        if (std::fabs(step_ - step_tmp) > 1e-7) {
          std::cerr << "Error while load_formatted of DiscreteAxis " << filename
                    << "\n";
          std::cerr << "step_ size is not consistent in line " << lines << "\n";
          std::cerr << value_tmp_new << " " << value_tmp_old << "\n";
          std::cerr << step_ << " vs. " << step_tmp << "\n";
        }
      }
      value_tmp_old = value_tmp_new;
      ++lines;
    }
  }
}

inline void serialize(const DiscreteAxis& t, std::ostream& out){
    serialize(t.min_, out);
    serialize(t.step_, out);
    serialize(t.invstep_, out);
    serialize(t.numBins_, out);
    serialize(t.max_, out);
    serialize(t.array_, out);
}

inline void deserialize(DiscreteAxis& t, std::istream& in){
    deserialize(t.min_, in);
    deserialize(t.step_, in);
    deserialize(t.invstep_, in);
    deserialize(t.numBins_, in);
    deserialize(t.max_, in);
    deserialize(t.array_, in);
}

} /* end of namespace betamc */ 
