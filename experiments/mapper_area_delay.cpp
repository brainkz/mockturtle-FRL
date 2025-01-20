/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/rsfq/rsfq_network_conversion.hpp>
#include <mockturtle/algorithms/rsfq/rsfq_path_balancing.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/algorithms/retiming.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/xmg.hpp>
#include <mockturtle/generators/arithmetic.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/rsfq_view.hpp>

#include <experiments.hpp>

constexpr int fDFF   = 0;
// constexpr int fNOT   = 1;
// constexpr int fMERGE = 2;
// constexpr int fOR    = 3;
// constexpr int fAND   = 4;
// constexpr int fXOR   = 5;
// constexpr int fOR3   = 6;
// constexpr int fAND3  = 7;
// constexpr int fMAJ3  = 8;
// constexpr int fCB    = 9;
constexpr int fSPL   = 10;
// constexpr int fPI    = 11;
// constexpr int fNOFUNC= 99;

// constexpr std::array<int,12> COSTS = {7, 9, 8, 8, 8, 7, 11, 11, 11, 8, 7, 0};
// constexpr std::array<int,12> COSTS_CONNECT = {6, 10, 7, 7, 7, 11, 999, 999, 999, 7, 3, 0};

// Removed input buffers in AND/OR gates
constexpr std::array<int,12> COSTS_CONNECT = {6, 9, 7, 3, 3, 11, 999, 999, 999, 7, 3, 0};

template<typename Ntk, typename AdderFn>
Ntk create_adder( uint32_t width, AdderFn&& adder )
{
  Ntk ntk;

  std::vector<typename Ntk::signal> a( width ), b( width );
  std::generate( a.begin(), a.end(), [&ntk]() { return ntk.create_pi(); } );
  std::generate( b.begin(), b.end(), [&ntk]() { return ntk.create_pi(); } );
  auto carry = ntk.get_constant( false );

  adder( ntk, a, b, carry );

  std::for_each( a.begin(), a.end(), [&]( auto f ) { ntk.create_po( f ); } );
  ntk.create_po( carry );

  assert( ntk.num_pis() == 2 * width );
  assert( ntk.num_pos() == width + 1 );

  return ntk;
}

typedef mockturtle::klut_network klut;
typedef mockturtle::xag_network xag;
typedef mockturtle::xmg_network xmg;
typedef mockturtle::aig_network aig;


template <typename Ntk>
std::tuple<mockturtle::binding_view<klut>, mockturtle::map_stats, int, int, bool> map_with_pb 
( 
  const std::string & benchmark, 
  const Ntk & tech_indep_ntk, 
  const mockturtle::tech_library<4u, mockturtle::classification_type::p_configurations> & tech_lib, 
  std::unordered_map<std::string, int> & nDFF_global, 
  bool area_oriented = false 
)
{
  mockturtle::map_params ps;
  ps.cut_enumeration_ps.minimize_truth_table = true;
  ps.cut_enumeration_ps.cut_limit = 24;
  ps.buffer_pis = false;
  if (area_oriented)
  {
      ps.skip_delay_round = true;
      ps.required_time = std::numeric_limits<float>::max();
  }
  mockturtle::map_stats st;
  mockturtle::binding_view<klut> res = map( tech_indep_ntk, tech_lib, ps, &st );
  mockturtle::depth_view<mockturtle::binding_view<klut>> dv { res };

  std::map<klut::node, int> dff_count;
  std::map<klut::node, int> fanout_count;

  /* RSFQ path balancing */
  auto balanced_res = mockturtle::rsfq_path_balancing( res );

  mockturtle::retime_params rps;
  mockturtle::retime_stats rst;
  auto net = mockturtle::rsfq_generic_network_create_from_mapped( balanced_res );
  mockturtle::retime( net, rps, &rst );
  auto retime_res = mockturtle::rsfq_mapped_create_from_generic_network( net );

  uint32_t num_ext_dffs = retime_res.num_dffs();
  
  uint32_t num_int_dffs = 0;
  retime_res.foreach_node( [&]( auto const& n ) {
    if ( !retime_res.has_binding( n ) )
      return;
    auto ndff = nDFF_global[retime_res.get_binding( n ).name];
    fmt::print("Adding {} internal DFFs from {}\n", ndff,  retime_res.get_binding( n ).name);
    num_int_dffs += ndff;
  } );

  /* RSFQ splitter insertion */
  uint32_t num_splitters = 0;
  retime_res.foreach_node( [&]( auto const& n ) {
    if ( !retime_res.is_constant( n ) )
      num_splitters += retime_res.fanout_size( n ) - 1;
  } );

  bool cec = rsfq_check_buffering( retime_res );
  cec &= benchmark == "hyp" ? true : experiments::abc_cec( retime_res, benchmark );

  // int internal_dff = 0;

  // dv.foreach_node( [&]( auto const& node ) 
  // {
  //     if ( dv.has_binding( node ) )
  //     {
  //       auto const& g = dv.get_binding( node );
        // fmt::print("Node {}\tGate {}\tnDFF {}\n", node, g.name, nDFF_global.at(g.name));
        // internal_dff += nDFF_global[g.name];//[];
        // dv.foreach_fanin(node, [&](auto pred, auto i) {
        //   int expected_lvl = dv.level(node) - g.pins[i].rise_block_delay;//g.num_vars-1-
        //   int actual_lvl = dv.level(pred);
        //   int extra_delay = expected_lvl - actual_lvl;

        //   dff_count[pred] = std::max(dff_count[pred], extra_delay);
          // fanout_count[pred]++;
        // });
  //     }
  // } );

  // int PB_DFF = std::accumulate(dff_count.begin(), dff_count.end(), 0, [](int a, const auto& p){ return a + p.second; });
  // PB_DFF += internal_dff;
  // int ext_SPL = std::accumulate(fanout_count.begin(), fanout_count.end(), 0, [](int a, const auto& p){ 
  //     return (p.second > 1)?(a + p.second - 1):(a); 
  //     });
  
  // Internal DFF area is already counted in the library
  int total_area = st.area + COSTS_CONNECT[fSPL] * num_splitters +  COSTS_CONNECT[fDFF] * num_ext_dffs;
  fmt::print("\t{} : Int: {}, Ext: {}, ratio: {}\n", benchmark, num_int_dffs, num_ext_dffs, (float)num_int_dffs / (num_int_dffs + num_ext_dffs) );
  return std::make_tuple( res, st, num_int_dffs + num_ext_dffs, total_area, cec );
}

// Function to read unordered_map from CSV file
std::unordered_map<std::string, int> readCSV(const std::string& filename) 
{
    std::ifstream infile(filename);             // Open the input file stream
    std::unordered_map<std::string, int> map;   // Create the unordered_map
    
    std::string line;
    std::getline(infile, line);                 // Ignore the header row

    fmt::print("READING CSV : {}\n", filename);
    // Read each subsequent row and add the key-value pair to the unordered_map
    while (std::getline(infile, line)) 
    {
        std::stringstream ss(line);
        std::string key;
        int value;
        std::getline(ss, key, ',');
        ss >> value;
        map[key] = value;
        fmt::print("READ {} : {}\n", key, value);
    }
    infile.close(); // Close the input file stream
    return map;
}

/* Gate costs are based on CONNECT library (from Japan) */
// const std::string DATABASE_PATH { "LIBRARY_2023_05_19_CONNECT.genlib" } ;
/* CONNECT library (from Japan) */
// const std::string DATABASE_PATH { "/Users/brainkz/Documents/GitHub/mockturtle_alessandro/build/LIBRARY_VANILLA_CONNECT.genlib" } ;
// LIBRARY_2023_05_19_CONNECT.genlib
// const std::string DATABASE_PATH { "/Users/brainkz/Documents/GitHub/mockturtle_alessandro/build/LIBRARY_2023_06_26_CONNECT_1111.genlib" } ;
// const std::string DATABASE_PATH { "/Users/brainkz/Documents/GitHub/mockturtle_alessandro/build/LIBRARY_2023_05_19_CONNECT.genlib" } ;
const std::string DATABASE_PATH { "/Users/brainkz/Documents/GitHub/mockturtle/build/LIBRARY_2023_06_27_CONNECT_CONSERVATIVE.genlib" } ;
/*The number of internal DFFs within each cell. 
Some of them are necessary not only for path balancing but also 
for synchronizing the pulses for AND gates. I include them 
in total DFF count */
// const std::string NDFF_PATH { "/Users/brainkz/Documents/GitHub/mockturtle_alessandro/build/nDFF_2023_05_08_CONNECT.csv" } ; 
// const std::string NDFF_PATH { "/Users/brainkz/Documents/GitHub/mockturtle/build/nDFF_2023_06_27_CONNECT_CONSERVATIVE.csv" } ; 
// const std::string NDFF_PATH { "/Users/brainkz/Documents/GitHub/mockturtle/build/NDFF_PARSED_2023_06_27_CONNECT_CONSERVATIVE.csv" } ; //updated costs here
const std::string NDFF_PATH { "/Users/brainkz/Documents/GitHub/mockturtle/build/NDFF_PARSED_CONNECT.csv" } ; //updated costs with filtering here


int main()
{
  using namespace experiments;
  using namespace mockturtle;
  
  
    std::unordered_map<std::string, std::tuple<double,double,double>> PBMAP;
    PBMAP["int2float"] = std::make_tuple(  270,   6432,  16);
    PBMAP["priority"]  = std::make_tuple( 9064, 102085, 127);
    PBMAP["sin"]       = std::make_tuple(13666, 215318, 182);
    PBMAP["cavlc"]     = std::make_tuple(  522,  16339,  17);
    PBMAP["dec"]       = std::make_tuple(    8,   5469,   4);
    PBMAP["c499"]      = std::make_tuple(  476,   7758,  13);
    PBMAP["c880"]      = std::make_tuple(  774,  12909,  22);
    PBMAP["c1908"]     = std::make_tuple(  696,  12013,  20);
    PBMAP["c3540"]     = std::make_tuple( 1159,  28300,  31);
    PBMAP["c5315"]     = std::make_tuple( 2908,  52033,  23);
    PBMAP["c7552"]     = std::make_tuple( 2429,  48482,  19);
    

  experiment<std::string, 
            double, double, double, 
            double, double, double, 
            double, double, double, 
            float> exp(
      "mapper", "benchmark", 
      "#DFF (base)", "#DFF (our)", "#DFF (ratio)", 
      "area (base)", "area (our)", "area (ratio)", 
      "delay (base)", "delay (our)", "delay (ratio)", 
      "time");

  fmt::print( "[i] processing technology library\n" );

  /* library to map to technology */
  std::vector<gate> gates;
  std::ifstream inputFile( DATABASE_PATH );
  std::unordered_map<std::string, int> nDFF_global = readCSV( NDFF_PATH );
  // std::unordered_map<std::string, int> nDFF_global;

  // Prints # DFF for each cell in the library
  // /*
    int ctr = 0;
    for (auto & [key, value]: nDFF_global)
    {
      fmt::print("{}:\t{}\t", key, value);
      if (++ctr % 8 == 0) 
      {
        std::cout << std::endl;
      }
    }
  // */
  
  if ( lorina::read_genlib( inputFile, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  mockturtle::tech_library_params tps; // tps.verbose = true;
  tech_library<4, mockturtle::classification_type::p_configurations> tech_lib( gates, tps );

  // // Where to find the MCNC benchmarks?
  auto benchmarks1 = epfl_benchmarks( experiments::adder | experiments::sin | experiments::cavlc | experiments::int2float | experiments::priority | experiments::i2c | experiments::voter | experiments::dec );
  // //   auto benchmarks1 = epfl_benchmarks( experiments::epfl & ~experiments::div & ~experiments::hyp & ~experiments::log2 & ~experiments::sqrt );
  auto benchmarks2 = iscas_benchmarks( experiments::c432 | experiments::c499 | experiments::c880 | experiments::c1355 | experiments::c1908 | experiments::c3540 | experiments::c5315 | experiments::c7552 );
  benchmarks1.insert(benchmarks1.end(), benchmarks2.begin(), benchmarks2.end());

  for ( auto const& benchmark : benchmarks1 )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    mockturtle::xmg_network tech_indep_ntk;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( tech_indep_ntk ) ) != lorina::return_code::success )
    {
      continue;
    }


    auto [res_area , st_area , PB_DFF_area , nJJ_area , cec_area ] = map_with_pb(benchmark, tech_indep_ntk, tech_lib, nDFF_global, true );
    auto [res_delay, st_delay, PB_DFF_delay, nJJ_delay, cec_delay] = map_with_pb(benchmark, tech_indep_ntk, tech_lib, nDFF_global, false);

    if (PBMAP.find(benchmark) != PBMAP.end())
    {
        const double ndff  = std::get<0>(PBMAP[benchmark]);
        const double njj   = std::get<1>(PBMAP[benchmark]);
        const double delay = std::get<2>(PBMAP[benchmark]);
        
//         exp( benchmark + " area" , ndff, PB_DFF_area,     PB_DFF_area/ndff, 
//                                    njj,   nJJ_area ,        nJJ_area/njj,
//                                    delay, st_area.delay, st_area.delay/delay,
//                                    to_seconds( st_area.time_total  )  );
        
        exp( benchmark + " delay", ndff, PB_DFF_delay,     PB_DFF_delay/ndff, 
                                   njj,   nJJ_delay ,        nJJ_delay/njj,
                                   delay, st_delay.delay, st_delay.delay/delay,
                                   to_seconds( st_delay.time_total * 1000 )  );
    }    
    else
    {        
//         exp( benchmark + " area" , 0xFFFFFFFF, PB_DFF_area,    0xFFFFFFFF, 
//                                    0xFFFFFFFF,   nJJ_area ,    0xFFFFFFFF,
//                                    0xFFFFFFFF, st_area.delay,  0xFFFFFFFF,
//                                    to_seconds( st_area.time_total  )  );
        
        exp( benchmark + " delay", 0xFFFFFFFF, PB_DFF_delay,   0xFFFFFFFF, 
                                   0xFFFFFFFF,   nJJ_delay ,   0xFFFFFFFF,
                                   0xFFFFFFFF, st_delay.delay, 0xFFFFFFFF,  
                                   to_seconds( st_delay.time_total * 1000 ) ); //
    }    
    // Write .bench file for inspection
    // mockturtle::write_bench(res_area , benchmark + "_connect_area"  + ".bench");
    mockturtle::write_bench(res_delay, benchmark + "_connect_full" + ".bench");
    exp.save();
    exp.table(); // print continuously
  }

  return 0;
}