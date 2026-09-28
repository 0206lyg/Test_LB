#ifndef SLURRY_GR_RE2_HISTORY_CSV_H
#define SLURRY_GR_RE2_HISTORY_CSV_H
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace slurry { namespace gr_re2 {
// Existing strict-solver histories contain no passed maximum-iteration steps.
// Migrate only the copied restart history, before opening it for append. No
// previously recorded physical quantity or sampling interval is changed.
inline void upgradePassMaxHistory(const std::filesystem::path& path) {
  std::ifstream source(path);
  std::string header;
  if(!std::getline(source,header))throw std::runtime_error("Cannot read restart history header");
  if(!header.empty()&&header.back()=='\r')header.pop_back();
  const std::string suffix=",max_iteration_passes,max_iteration_passes_total";
  if(header.size()>=suffix.size()&&header.compare(header.size()-suffix.size(),suffix.size(),suffix)==0)return;
  if(header.find("max_iteration_passes")!=std::string::npos)
    throw std::runtime_error("Incompatible maximum-iteration history columns");
  const auto temporary=std::filesystem::path(path.string()+".pass-max.tmp");
  try {
    std::ofstream target(temporary,std::ios::trunc);
    if(!target)throw std::runtime_error("Cannot migrate restart history");
    target<<header<<suffix<<'\n';
    std::string row;
    while(std::getline(source,row)){
      if(!row.empty()&&row.back()=='\r')row.pop_back();
      if(!row.empty())target<<row<<",0,0\n";
    }
    if(source.bad())throw std::runtime_error("Cannot finish reading restart history");
    target.close();source.close();
    if(!target)throw std::runtime_error("Cannot finish migrating restart history");
    std::filesystem::rename(temporary,path);
  } catch(...) {
    std::error_code ignored;std::filesystem::remove(temporary,ignored);throw;
  }
}
} }
#endif
