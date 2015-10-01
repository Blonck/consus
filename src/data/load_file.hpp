#pragma once

#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_parse.hpp>
#include <boost/spirit/include/qi_numeric.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/progress.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream.hpp>

#include "../vector/helper.hpp"

namespace consus
{

namespace impl {

bool check_and_sort_columns(std::vector<size_t>& use_cols, size_t max_size) {
  auto it = std::find_if(
      use_cols.begin(), use_cols.end(),
      [&](const size_t& elem) -> bool { return elem >= max_size; });
  if (it != use_cols.end()) {
    std::cerr << "Column " << *it << " is not a valid columns\n";
    return false;
  }
  std::sort(use_cols.begin(), use_cols.end());
  if (std::adjacent_find(use_cols.begin(), use_cols.end()) != use_cols.end()) {
    std::cerr << "Given two identical columns!\n";
    return false;
  }
  return true;
}

}

// TODO: refactor read_ssv and read_ssv_jk, nearly identical -> identical
//       parts should use same functions
// TODO: read_header is not necessary, read_ssv should return std::pair<vec2d, header>
// TODO: if the first line which is not a comment is parsed by hand, one can estimate
//       the number of elements and do a reserve on the corresponding vector
// TODO: atm only ascii and .gz data can be read in, no problem to add bz2 or something else

/// read file with space separated values
///
/// all line which starts with # will be ignored
/// data will be transposed so that the return value is a
/// 2 dim vector, where every vector contains 1 row, by convention each row is a
/// single time series
/// this function allocates temporarily two times the size of the file as memory
/// \param filename - name of the file which will be read in and parsed
/// TODO: parse atm every value as double and convert to T: should be changed
template <class T = double>
std::vector<std::vector<T>> read_ssv(const std::string& filename,
                                     std::vector<size_t> use_cols = {}) {
  std::string content;
  if (boost::iends_with(filename, ".gz")) {
    std::ifstream infile(filename, std::ios::in | std::ios::binary);
    if (infile) {
      boost::iostreams::filtering_streambuf<boost::iostreams::input> buf;
      buf.push(boost::iostreams::gzip_decompressor());
      buf.push(infile);
      std::stringstream ss(std::ios::in | std::ios::out | std::ios::binary);
      boost::iostreams::copy(buf, ss);
      content = ss.str();
    } else {
      std::cerr << filename << "does not exists\n";
      std::exit(1);
    }
  } else {
    std::ifstream infile(filename, std::ios::in);
    if (infile) {
      infile.seekg(0, std::ios::end);
      content.resize(infile.tellg());
      infile.seekg(0, std::ios::beg);
      infile.read(&content[0], content.size());
      infile.close();
    } else {
      std::cerr << filename << "does not exists\n";
      std::exit(1);
    }
  }

  using boost::spirit::qi::double_;
  using boost::spirit::qi::eol;
  using boost::spirit::qi::repeat;
  using boost::spirit::qi::phrase_parse;
  using boost::spirit::qi::lexeme;
  using boost::spirit::ascii::space;
  using boost::spirit::ascii::blank;
  using boost::spirit::ascii::char_;
  std::string::const_iterator startit = content.begin();
  std::string::const_iterator endit = content.end();
  std::vector<std::vector<T>> v;
  bool r = phrase_parse(startit, endit, (repeat[(repeat[double_] >> eol)]),
                        blank | '#' >> *(char_ - eol) >> eol, v);
  if(not r or startit != endit){
    std::cerr << "error on parsing file" << filename << "\n";
    std::exit(1);
  }
  size_t columns = v[0].size();
  for(size_t i=0; i < v.size(); ++i){
    if (columns != v[i].size()){
      std::cerr << "error in file " << filename << " on line " << i << "\n";
      std::exit(1);
    }
  }
  content.clear();

  std::vector<std::vector<T>> vT;
  // transpose data
  if (use_cols.empty()) {
    vT.resize(v[0].size(), vec1d(v.size()));
    for (size_t i = 0; i < v.size(); ++i) {
      for (size_t j = 0; j < v[j].size(); ++j) {
        vT[j][i] = v[i][j];
      }
    }
  }else{
    // if not all columns are valid stop the program
    if (not impl::check_and_sort_columns(use_cols, v[0].size())){
      std::exit(1);
    }
    vT.resize(use_cols.size(), vec1d(v.size()));
    for (size_t i = 0; i < v.size(); ++i) {
      for (size_t j = 0; j < use_cols.size(); ++j) {
        vT[j][i] = v[i][use_cols[j]];
      }
    }
  }
  return vT;
}

template <class T = double>
std::vector<std::vector<T>> read_ssv_jk(const std::string& filename,
                                        size_t num_bin, size_t num_bins,
                                        std::vector<size_t> use_cols = {}) {
  std::string content;
  if (boost::iends_with(filename, ".gz")) {
    std::ifstream infile(filename, std::ios::in | std::ios::binary);
    if (infile) {
      boost::iostreams::filtering_streambuf<boost::iostreams::input> buf;
      buf.push(boost::iostreams::gzip_decompressor());
      buf.push(infile);
      std::stringstream ss(std::ios::in | std::ios::out | std::ios::binary);
      boost::iostreams::copy(buf, ss);
      content = ss.str();
    } else {
      std::cerr << filename << "does not exists\n";
      std::exit(1);
    }
  } else {
    std::ifstream infile(filename, std::ios::in);
    if (infile) {
      infile.seekg(0, std::ios::end);
      content.resize(infile.tellg());
      infile.seekg(0, std::ios::beg);
      infile.read(&content[0], content.size());
      infile.close();
    } else {
      std::cerr << filename << "does not exists\n";
      std::exit(1);
    }
  }

  using boost::spirit::qi::double_;
  using boost::spirit::qi::eol;
  using boost::spirit::qi::repeat;
  using boost::spirit::qi::phrase_parse;
  using boost::spirit::qi::lexeme;
  using boost::spirit::ascii::space;
  using boost::spirit::ascii::blank;
  using boost::spirit::ascii::char_;
  std::string::const_iterator startit = content.begin();
  std::string::const_iterator endit = content.end();
  std::vector<std::vector<T>> v;
  bool r = phrase_parse(startit, endit, (repeat[(repeat[double_] >> eol)]),
                        blank | '#' >> *(char_ - eol) >> eol, v);
  if(not r or startit != endit){
    std::cerr << "error on parsing file" << "\n";
    std::abort();
  }
  size_t columns = v[0].size();
  for(size_t i=0; i < v.size(); ++i){
    if (columns != v[i].size()){
      std::cerr << "error on line " << i << "\n";
      std::abort();
    }
  }
  content.clear();
  std::vector<std::vector<T>> v_jk;
  size_t length_bins = v.size() / num_bins;
  if (v.size() % num_bins != 0) {
    std::cerr << "bin size and time series size doesn't fit\n";
    std::abort();
  }
  v_jk.reserve(v.size()-length_bins);
  for (size_t k = 0; k < num_bin; ++k) {
    for (size_t l = 0; l < length_bins; ++l) {
      v_jk.push_back(v[k * length_bins + l]);
    }
  }
  for (size_t k = num_bin + 1; k < num_bins; ++k) {
    for (size_t l = 0; l < length_bins; ++l) {
      v_jk.push_back(v[k * length_bins + l]);
    }
  }
  v_jk.shrink_to_fit();
  v = std::move(v_jk);

  std::vector<std::vector<T>> vT;
  if (use_cols.empty()) {
    vT.resize(v[0].size(), vec1d(v.size()));
    for (size_t i = 0; i < v.size(); ++i) {
      for (size_t j = 0; j < v[j].size(); ++j) {
        vT[j][i] = v[i][j];
      }
    }
    return vT;
  }else{
    if (not impl::check_and_sort_columns(use_cols, v[0].size())){
      std::exit(1);
    }
    vT.resize(use_cols.size(), vec1d(v.size()));
    for (size_t i = 0; i < v.size(); ++i) {
      for (size_t j = 0; j < use_cols.size(); ++j) {
        vT[j][i] = v[i][use_cols[j]];
      }
    }
    return vT;
  }
}

/// parse the header of a file containing time series
///
/// returns the space separated entries of the last line which is commented out
/// with # before the time series starts
/// return value - std::vector<std::string>
/// \param filename - name of the file which is be parsed
std::vector<std::string> read_header(const std::string& filename,
                                     std::vector<int> use_cols = {}) {
  std::string content;
  std::stringstream data;
  if (boost::iends_with(filename, ".gz")) {
    std::ifstream infile(filename, std::ios::in | std::ios::binary);
    if (infile) {
      boost::iostreams::filtering_streambuf<boost::iostreams::input> buf;
      buf.push(boost::iostreams::gzip_decompressor());
      buf.push(infile);
      std::stringstream ss(std::ios::in | std::ios::out | std::ios::binary);
      boost::iostreams::copy(buf, ss);
      content = ss.str();
    } else {
      std::cerr << filename << "does not exists\n";
      std::exit(1);
    }
  } else {
    std::ifstream infile(filename, std::ios::in);
    if (infile) {
      infile.seekg(0, std::ios::end);
      content.resize(infile.tellg());
      infile.seekg(0, std::ios::beg);
      infile.read(&content[0], content.size());
      infile.close();
    } else {
      std::cerr << filename << "does not exists\n";
      std::exit(1);
    }
  }
  std::string line;
  std::string lastline;
  data << content;
  while( std::getline(data, line) ){
    if ( line[0] == '#'){
      lastline = line;
    }else{
      break;
    }
  }

  using boost::spirit::qi::double_;
  using boost::spirit::qi::eol;
  using boost::spirit::qi::repeat;
  using boost::spirit::qi::phrase_parse;
  using boost::spirit::qi::lexeme;
  using boost::spirit::ascii::space;
  using boost::spirit::ascii::blank;
  using boost::spirit::ascii::char_;
  std::string::const_iterator startit = lastline.begin() + 1;
  std::string::const_iterator endit = lastline.end();
  std::vector<std::string> header;
  auto r = phrase_parse(startit, endit, (repeat[lexeme[+(char_ - blank)]]), blank,
                   header);
  if(not r or startit != endit){
    std::cerr << "error" << "\n";
    std::cout << std::endl;
    std::abort();
  }
  if( use_cols.size() != 0 ){
    std::vector<std::string> new_header;
    new_header.reserve(use_cols.size());
    for(size_t i = 0; i < use_cols.size(); ++i){
      new_header.push_back( header[use_cols[i]] );
    }
    return new_header;
  }
  return header;
}

} /* end of namespace consus */ 
