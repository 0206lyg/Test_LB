/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef SLURRY_GR_RE2_COUPLED_CHECKPOINT_H
#define SLURRY_GR_RE2_COUPLED_CHECKPOINT_H
#include "quickConfig.h"
#include "particleSubsteps.h"
#include <array>
#include <cstring>
#include <iomanip>
#include <map>
#include <sstream>
#include <type_traits>

// SLURRY SCOPE BEGIN
namespace slurry { namespace gr_re2 { namespace checkpoint_detail {
namespace fs=std::filesystem;
using U64=std::uint64_t;
constexpr U64 format=1,byteOrder=0x0102030405060708ULL;
constexpr U64 hashBasis=14695981039346656037ULL;
inline void hashBytes(U64& hash,const char* data,std::size_t size) {
  for(std::size_t i=0;i<size;++i){hash^=static_cast<unsigned char>(data[i]);hash*=1099511628211ULL;}
}
class Output {
  std::ofstream stream;U64 hash=hashBasis;
public:
  explicit Output(const fs::path& path):stream(path,std::ios::binary) {
    if(!stream)throw std::runtime_error("Cannot create checkpoint: "+path.string());
  }
  void bytes(const void* data,std::size_t size) {
    stream.write(static_cast<const char*>(data),size);hashBytes(hash,static_cast<const char*>(data),size);
    if(!stream)throw std::runtime_error("Checkpoint write failed");
  }
  template<class V>void put(const V& value) {
    static_assert(std::is_trivially_copyable<V>::value,"Checkpoint values must not contain pointers");
    bytes(&value,sizeof(V));
  }
  void finish(){stream.write(reinterpret_cast<const char*>(&hash),sizeof(hash));stream.close();
    if(!stream)throw std::runtime_error("Checkpoint close failed");}
};
// Verify the complete file before allowing a loader to modify any live state.
class Input {
  std::ifstream stream;U64 remaining=0;
public:
  explicit Input(const fs::path& path):stream(path,std::ios::binary) {
    if(!stream)throw std::runtime_error("Cannot read checkpoint: "+path.string());
    const auto size=fs::file_size(path);
    if(size<sizeof(U64))throw std::runtime_error("Truncated checkpoint: "+path.string());
    remaining=size-sizeof(U64);U64 left=remaining,hash=hashBasis,expected=0;
    std::array<char,65536> buffer{};
    while(left){const auto count=std::min<U64>(left,buffer.size());stream.read(buffer.data(),count);
      if(!stream)throw std::runtime_error("Truncated checkpoint: "+path.string());
      hashBytes(hash,buffer.data(),count);left-=count;}
    stream.read(reinterpret_cast<char*>(&expected),sizeof(expected));
    if(!stream||hash!=expected)throw std::runtime_error("Checkpoint checksum mismatch: "+path.string());
    stream.clear();stream.seekg(0);
  }
  void bytes(void* data,std::size_t size) {
    if(size>remaining)throw std::runtime_error("Checkpoint payload length mismatch");
    stream.read(static_cast<char*>(data),size);remaining-=size;
    if(!stream)throw std::runtime_error("Truncated checkpoint payload");
  }
  template<class V>void get(V& value){static_assert(std::is_trivially_copyable<V>::value,"");bytes(&value,sizeof(V));}
  template<class V>void expect(const V& expected){V actual{};get(actual);
    if(actual!=expected)throw std::runtime_error("Checkpoint format/ABI mismatch; use the same OpenLB layout");}
  void finish(){if(remaining)throw std::runtime_error("Unexpected trailing checkpoint payload");}
};
inline std::map<std::string,double> physics(const Config& c,const Units& u,std::size_t count,int ranks) {
  std::map<std::string,double> out;
#define KEEP(k) out[#k]=c.k
  KEEP(shear_rate);KEEP(box_x);KEEP(box_y);KEEP(box_z);KEEP(dx);KEEP(diameter);KEEP(thickness);
  KEEP(rho_particle);KEEP(rho_fluid);KEEP(dynamic_viscosity);KEEP(nu_lattice);KEEP(epsilon_cells);
  KEEP(hamaker);KEEP(sigma_lj);KEEP(switch_gap);KEEP(cutoff_gap);
  if(c.surface_adhesion){
    KEEP(surface_adhesion);KEEP(adhesion_work);KEEP(adhesion_range);
    KEEP(curvature_switch_gap);KEEP(curvature_cutoff_gap);
    out["interaction_model_version"]=1.;
  } else {
    // Preserve the old signature exactly so explicit legacy restarts stay valid.
    KEEP(local_gap);KEEP(local_gap_fraction);KEEP(local_switch_excess_gap);KEEP(local_cutoff_excess_gap);
  }
  KEEP(lubrication_cutoff_cells);
  KEEP(rough_contact_enabled);KEEP(roughness_gap);KEEP(sliding_friction);KEEP(tangential_stiffness);
  KEEP(rolling_length);KEEP(rolling_yield_angle);
#undef KEEP
  out["dt_s"]=u.dt;out["particle_count"]=count;out["ranks"]=ranks;
  return out;
}
inline std::string signature(const Config& c,const Units& u,std::size_t count,int ranks) {
  std::ostringstream out;out<<std::setprecision(17);
  for(const auto& entry:physics(c,u,count,ranks))out<<entry.first<<'='<<entry.second<<'\n';
  return out.str();
}
struct State {
  U64 step=0;
  std::vector<graphite::Body> bodies;
  std::vector<graphite::GapCache> cache;
  graphite::PersistentContactState contacts;
  std::vector<graphite::Vec3> angularAcceleration;
  graphite::ParticleStepDiagnostics diagnostic;
  // Accumulated wall, fluid, mapping, coupling, particle, output seconds.
  std::array<double,6> timing{};
};
inline void saveState(const fs::path& path,const State& s,const std::string& expected) {
  const auto count=s.bodies.size(),pairs=count*(count-1)/2;
  if(s.cache.size()!=pairs||s.contacts.size()!=pairs||s.angularAcceleration.size()!=count)
    throw std::runtime_error("Incomplete particle/contact checkpoint state");
  Output out(path);out.put(U64(0x4752535441544531ULL));out.put(format);out.put(byteOrder);
  for(U64 size:{sizeof(graphite::Body),sizeof(graphite::GapCache),sizeof(graphite::RoughContactState),
                sizeof(graphite::ParticleStepDiagnostics)})out.put(size);
  out.put(U64(expected.size()));out.bytes(expected.data(),expected.size());
  out.put(s.step);out.put(U64(s.bodies.size()));out.put(s.diagnostic);out.put(s.timing);
  for(const auto& body:s.bodies)out.put(body);
  for(const auto& cache:s.cache)out.put(cache);
  for(const auto& contact:s.contacts)out.put(contact);
  for(const auto& acceleration:s.angularAcceleration)out.put(acceleration);
  out.finish();
}
inline State loadState(const fs::path& path,std::size_t count,const std::string& expected) {
  Input in(path);in.expect(U64(0x4752535441544531ULL));in.expect(format);in.expect(byteOrder);
  for(U64 size:{sizeof(graphite::Body),sizeof(graphite::GapCache),sizeof(graphite::RoughContactState),
                sizeof(graphite::ParticleStepDiagnostics)})in.expect(size);
  U64 length=0;in.get(length);if(length!=expected.size())throw std::runtime_error("Checkpoint physical settings changed");
  std::string saved(length,'\0');in.bytes(&saved[0],length);
  if(saved!=expected)throw std::runtime_error("Checkpoint physical settings, timestep, particle count or MPI ranks changed");
  State s;in.get(s.step);in.expect(U64(count));in.get(s.diagnostic);in.get(s.timing);
  const auto pairs=count*(count-1)/2;
  s.bodies.resize(count);s.cache.resize(pairs);s.contacts.resize(pairs);s.angularAcceleration.resize(count);
  for(auto& body:s.bodies)in.get(body);
  for(auto& cache:s.cache)in.get(cache);
  for(auto& contact:s.contacts)in.get(contact);
  for(auto& acceleration:s.angularAcceleration)in.get(acceleration);
  in.finish();return s;
}
template<class Lattice>void saveLattice(Lattice& lattice,const fs::path& path) {
  Output out(path);out.put(U64(0x47524c4154544931ULL));out.put(format);out.put(byteOrder);
  olb::Serializer serial(lattice);serial.computeSize();out.put(U64(serial.getSize()));
  std::size_t size=0;U64 written=0;
  while(const bool* block=serial.getNextBlock(size,false)){out.bytes(block,size);written+=size;}
  if(written!=serial.getSize())throw std::runtime_error("Lattice checkpoint serialization size mismatch");
  out.finish();
}
template<class Lattice>void loadLattice(Lattice& lattice,const fs::path& path) {
  Input in(path);in.expect(U64(0x47524c4154544931ULL));in.expect(format);in.expect(byteOrder);
  olb::Serializer serial(lattice);serial.computeSize();in.expect(U64(serial.getSize()));
  std::size_t size=0;
  while(bool* block=serial.getNextBlock(size,true))in.bytes(block,size);
  in.finish();lattice.postLoad();
}
inline std::string rankFile(int rank){return "lattice_rank_"+std::to_string(rank)+".bin";}
inline std::string stepName(U64 step){std::ostringstream s;s<<"checkpoint_"<<std::setw(20)<<std::setfill('0')<<step;return s.str();}

template<class Lattice>fs::path save(Lattice& lattice,const Config& c,const Units& u,const State& state) {
  const int rank=olb::singleton::mpi().getRank(),ranks=olb::singleton::mpi().getSize();
  const fs::path root=fs::path(c.output_dir)/"checkpoints",final=root/stepName(state.step);
  const fs::path staging=final.string()+".partial";
  if(rank==0){fs::create_directories(root);
    if(fs::exists(final))throw std::runtime_error("Checkpoint step already exists: "+final.string());
    if(fs::exists(staging))fs::remove_all(staging);
    fs::create_directory(staging);}
  olb::singleton::mpi().barrier();
  saveLattice(lattice,staging/rankFile(rank));
  if(rank==0){saveState(staging/"state.bin",state,signature(c,u,state.bodies.size(),ranks));
    fs::copy_file(c.particles_csv,staging/"initial_particles.csv");}
  // Publish only after every rank has closed its lattice file successfully.
  olb::singleton::mpi().barrier();
  if(rank==0){
    std::vector<std::string> files{"state.bin","initial_particles.csv"};
    for(int i=0;i<ranks;++i)files.push_back(rankFile(i));
    std::ofstream meta(staging/"checkpoint.json");meta<<std::setprecision(17);
    meta<<"{\"format_version\":1,\"engine\":\"pure_gr\",\"complete\":true,\"step\":"<<state.step
        <<",\"time_s\":"<<state.step*u.dt<<",\"strain\":"<<state.step*u.dt*c.shear_rate
        <<",\"dt_s\":"<<u.dt<<",\"shear_rate_s_inv\":"<<c.shear_rate<<",\"ranks\":"<<ranks
        <<",\"immutable_config\":{";
    bool first=true;for(const auto& entry:physics(c,u,state.bodies.size(),ranks)){
      if(!first)meta<<',';
      first=false;meta<<std::quoted(entry.first)<<':'<<entry.second;}
    meta<<"},\"files\":{";first=true;
    for(const auto& file:files){if(!first)meta<<',';first=false;
      meta<<std::quoted(file)<<':'<<fs::file_size(staging/file);}
    meta<<"},\"output_bytes\":{";first=true;
    for(const std::string file:{"history.csv","particles.csv"}){
      if(!first)meta<<',';
      first=false;meta<<std::quoted(file)<<':'
        <<(fs::exists(fs::path(c.output_dir)/file)?fs::file_size(fs::path(c.output_dir)/file):0);}
    meta<<"}}\n";meta.close();if(!meta)throw std::runtime_error("Checkpoint metadata write failed");
    fs::rename(staging,final);
    const auto pointer=fs::path(c.output_dir)/"latest_checkpoint.txt";
    std::ofstream latest(pointer.string()+".tmp");latest<<"checkpoints/"<<final.filename().string()<<'\n';latest.close();
    if(!latest)throw std::runtime_error("Checkpoint pointer write failed");
    fs::rename(pointer.string()+".tmp",pointer);
    std::vector<fs::path> completed;
    for(const auto& entry:fs::directory_iterator(root)){
      const auto name=entry.path().filename().string();
      if(entry.is_directory()&&name.size()==31&&name.compare(0,11,"checkpoint_")==0
          &&name.find_first_not_of("0123456789",11)==std::string::npos
          &&fs::is_regular_file(entry.path()/"checkpoint.json"))completed.push_back(entry.path());}
    std::sort(completed.begin(),completed.end());
    while(completed.size()>c.checkpoint_keep){fs::remove_all(completed.front());completed.erase(completed.begin());}
  }
  olb::singleton::mpi().barrier();return final;
}
} } } // SLURRY SCOPE END
#endif
