#include "slurry/gr_re2/historyCsv.h"
#include <cassert>
#include <iostream>
#include <sstream>
using namespace slurry::gr_re2;
std::string read(const std::filesystem::path&p){std::ifstream in(p);std::ostringstream out;out<<in.rdbuf();return out.str();}
int main(){
 const auto p=std::filesystem::temp_directory_path()/"gr_history_pass_max_test.csv";
 {std::ofstream f(p);f<<"step,stress\n0,0\n20,12.5\n";}
 upgradePassMaxHistory(p);
 const auto result=read(p);
 assert(result=="step,stress,max_iteration_passes,max_iteration_passes_total\n0,0,0,0\n20,12.5,0,0\n");
 upgradePassMaxHistory(p);assert(read(p)==result);
 {std::ofstream f(p);f<<"step,max_iteration_passes\n0,0\n";}
 bool rejected=false;try{upgradePassMaxHistory(p);}catch(const std::runtime_error&){rejected=true;}
 assert(rejected);std::filesystem::remove(p);
 std::cout<<"PASS history migration, idempotence, incompatible-header rejection\n";
}
