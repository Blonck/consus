#pragma once

#include <vector>
#include <string>
#include "boost/filesystem.hpp"
#include "boost/regex.hpp"

std::vector<std::string> all_matching_paths(const std::string& target_path,
                                            const boost::regex& filter) {
  std::vector<std::string> matching_paths;

  boost::filesystem::directory_iterator end_itr;
  for ( boost::filesystem::directory_iterator i( target_path ); i != end_itr; ++i){
    if (not boost::filesystem::is_directory( i->status())){
      continue;
    }
    boost::smatch what;
    if (not boost::regex_match( i->path().string(), what, filter)){
      continue;
    }
    matching_paths.push_back(i->path().string());
  }
  return matching_paths;
}
