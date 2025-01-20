#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <tuple>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <stack>
#include <array>
#include <fmt/format.h>
#include <thread>
#include <kitty/kitty.hpp>
#include <parallel_hashmap/phmap.h>

typedef unsigned short US;
typedef unsigned int   UI;
typedef unsigned long  UL;
typedef const unsigned short CUS;
typedef const unsigned int   CUI;
typedef const unsigned long  CUL;
typedef unsigned long long ULL;

bool DEBUG = false;

constexpr UI NUM_VARS = 4u;
typedef kitty::static_truth_table<4u> TT4;
typedef kitty::static_truth_table<3u> TT3;
typedef kitty::static_truth_table<2u> TT2;
typedef kitty::static_truth_table<1u> TT1;
// typedef kitty::dynamic_truth_table    DTT;

constexpr US fDFF   = 0;
constexpr US fNOT   = 1;
constexpr US fMERGE = 2;
constexpr US fOR  = 3;
constexpr US fAND   = 4;
constexpr US fXOR   = 5;
constexpr US fOR3   = 6;
constexpr US fAND3  = 7;
constexpr US fMAJ3  = 8;
constexpr US fCB  = 9;
constexpr US fSPL   = 10;
constexpr US fPI  = 11;
constexpr US fT1  = 12;
constexpr US fNOFUNC= 99;

// constexpr std::array<int,12> COSTS_CONNECT = {6, 10, 7, 3, 3, 11, 999, 999, 999, 7, 3, 0};
// constexpr std::array<int,12> COSTS_CONNECT_UPD = {5, 9, 7, 3, 3, 11, 999, 999, 999, 7, 3, 0};
// constexpr std::array<int,12> COSTS_CONNECT_CONSERVATIVE = {6, 9, 7, 3, 3, 11, 999, 999, 999, 7, 3, 0};
constexpr std::array<int,12> COSTS_CONNECT_CONSERVATIVE = {6, 9, 7, 3, 3, 11, 999, 999, 999, 7, 3, 0};
constexpr std::array<int,12> COSTS_CONNECT_PESSIMISTIC_3IN = {6, 9, 7, 7, 7, 11, 10, 10, 10, 7, 3, 0};
constexpr std::array<int,12> COSTS_SUNMAGNETICS = {7, 9, 8, 8, 12, 8, 999, 999, 999, 8, 3, 0};
constexpr std::array<int,13> COSTS_SUNMAGNETICS_EXTENDED = {7, 9, 8, 8, 12, 8, 999, 999, 999, 8, 3, 0, 25};

constexpr UI kNumThreads = 100;

const std::string GENLIB_PHASE = "UNKNOWN";
constexpr float GENLIB_INPUT_LOAD = 1;
constexpr float GENLIB_MAX_LOAD = 999;
constexpr float GENLIB_RISE_BLOCK_DELAY   = 0.025;
constexpr float GENLIB_RISE_FANOUT_DELAY  = 0.025;
constexpr float GENLIB_FALL_BLOCK_DELAY   = 0.025;
constexpr float GENLIB_FALL_FANOUT_DELAY  = 0.025;


// constexpr bool accel_cost = true;
#define accel_cost false

#define VECTOR_CONTAINS(vec, value) (std::find(vec.begin(), vec.end(), value) != vec.end())

constexpr std::array<UI,12> COSTS = {7, 9, 8, 8, 8, 7, 11, 11, 11, 8, 7, 0}; // ORIGINAL COSTS
// constexpr std::array<UI,12> COSTS = {5,  9, 7, 3, 3, 11, 999, 999, 999, 8, 7, 0};
                  // {6, 10, 7, 7, 7, 11, 999, 999, 999, 7, 3, 0};
phmap::flat_hash_map<US, std::string> F2STR { 
  {fDFF   , "DFF"},
  {fNOT   , "NOT"},
  {fMERGE , "MRG"},
  {fOR    , "OR "},
  {fAND   , "AND"},
  {fXOR   , "XOR"},
  {fOR3   , "OR3"},
  {fAND3  , "AND3"},
  {fMAJ3  , "MAJ3"},
  {fCB    , "CB "},
  {fSPL   , "SPL"},
  {fPI    , "PI "},
  {fNOFUNC, "N/A"}
  }; 

std::array<uint16_t, NUM_VARS> PI_WORDS = {{0x5555, 0x3333, 0x0F0F, 0x00FF}};

const phmap::flat_hash_map<uint16_t, uint8_t> PIFUNC2IDX = {{0x5555, 0}, {0x3333, 1}, {0x0F0F, 2}, {0x00FF, 3}};

phmap::flat_hash_map<uint16_t, std::string> PI2LETTER { 
  {0x00FF, "d"},
  {0x0F0F, "c"},
  {0x3333, "b"},
  {0x5555, "a"},
  {0x0000, "0"},
  {0xFFFF, "1"},
  }; 

phmap::flat_hash_map<uint16_t, std::string> IDX2LETTER { 
  {3, "d"},
  {2, "c"},
  {1, "b"},
  {0, "a"},
  }; 

phmap::flat_hash_map<uint16_t, uint16_t> dummy_map {
  {0x00FF, 0x00FF},
  {0x0F0F, 0x0F0F},
  {0x3333, 0x3333},
  {0x5555, 0x5555}
  };

phmap::flat_hash_map<uint16_t, std::string> FN2LETTER { 
  {0x00FF, "d"},
  {0x0F0F, "c"},
  {0x3333, "b"},
  {0x5555, "a"},
  {0x0F, "c"},
  {0x33, "b"},
  {0x55, "a"},
  {0x3, "b"},
  {0x5, "a"},
  {0x1, "a"},
  {0x0000, "CONST0"},
  {0xFFFF, "CONST1"},
  }; 

constexpr UL NUM_TT = (1 << (1 << NUM_VARS));
constexpr UI ONES = NUM_TT - 1;
// constexpr UL INF = 0xFFFFFF;
constexpr UL INF = 0xFFFF;

// Hash combiner
template <typename T>
static void hash_combine(std::size_t& seed, const T& val) {
  seed ^= std::hash<T>{}(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

ULL calculate_hash(UI func, US last_func, UI cost, UI depth, bool xorable, std::vector<ULL> parent_hashes = {})
{
  std::size_t seed = 0;
  hash_combine(seed, func);
  hash_combine(seed, last_func);
  hash_combine(seed, cost);
  hash_combine(seed, depth);
  hash_combine(seed, xorable);
  for (const auto parent_hash : parent_hashes) {
    hash_combine(seed, parent_hash);
  }
  return seed;
};

// TODO: find P representative and minimize functions in support
// TODO: permute the delays accordingly

// TODO: in next function, test whether one delay pattern dominates the other. 
// TODO: Perhaps, a delay object will help


template <typename TT>
bool _tt_gt (const std::pair<TT, uint8_t>& a, const std::pair<TT, uint8_t>& b) { return ~(a.first < b.first); }

class Node 
{
public:
  UI func = 0;
  US last_func = fNOFUNC;
  UI cost = INF;
  UI depth = INF;
  bool xorable = false;
  std::vector<ULL> parent_hashes;
  UI lvl = INF;
  ULL hash = 0;

  // UI func = 0;                    // 16
  // US last_func = fNOFUNC;         // 7
  // UI cost = INF;                  // 16
  // UI depth = INF;                 // 5
  // bool xorable = false;           // 1
  // UI lvl = INF;                   // 3
  // ULL hash = 0;                   // 64
  // std::vector<ULL> parent_hashes; // 64 * vector_size

  Node() = default;
  Node(const Node& other) : func(other.func), last_func(other.last_func), cost(other.cost), depth(other.depth), xorable(other.xorable), parent_hashes(other.parent_hashes), hash(other.hash), lvl(other.lvl) {}

  Node& operator=(const Node& other) {
    if (this != &other) {
      func = other.func;
      last_func = other.last_func;
      cost = other.cost;
      depth = other.depth;
      xorable = other.xorable;
      parent_hashes = other.parent_hashes;
      hash = other.hash;
      lvl = other.lvl;
    }
    return *this;
  }

  Node(UI _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes)
    : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes(_parent_hashes), lvl(depth / 3)
  {
    // Calculate hash based on the hashes of parent_hashes and specified fields
    hash = calculate_hash(func, last_func, cost, depth, xorable, parent_hashes);
    // hash = calculate_hash();
  }

  Node(UI _func, US _last_func, UI _cost, UI _depth, bool _xorable)
    : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes{}, lvl(depth / 3)
  {
    // Calculate hash based on the hashes of parent_hashes and specified fields
    // It is assumed the node has no parents
    hash = calculate_hash(func, last_func, cost, depth, xorable);
  }

  Node(UI _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes, ULL _hash)
    : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes(_parent_hashes), lvl(depth / 3), hash(_hash)
  {
    // Hash is precalculated based on the hashes of parent_hashes and specified fields
  }

  ~Node() { }

  // Equality operator
  // bool operator==(const Node& other) const 
  // {
  //   return hash == other.hash;
  // }
  bool operator==(const Node& other) const 
  {
    return  func          == other.func &&
            last_func     == other.last_func &&
            cost          == other.cost &&
            depth         == other.depth &&
            xorable       == other.xorable &&
            parent_hashes == other.parent_hashes &&
            lvl           == other.lvl &&
            hash          == other.hash;
  }

  bool operator!=(const Node& other) const 
  {
    return !(*this == other);
  }

  ULL get_hash() 
  {
    return hash;
  }

  char pi_letter() 
  {
    if      (func == 0x5555) { return 'a'; }
    else if (func == 0x3333) { return 'b'; }
    else if (func == 0x0F0F) { return 'c'; }
    else if (func == 0x00FF) { return 'd'; }
    else 
    {
      throw fmt::format("The node is not a PI: {0}", to_str());
    }
  }

  std::string to_str() const
  {
    // std::string s;
    // for (auto & ph : parent_hashes)
    // {
    //   Node & ref = GNM[ph];
    //   s += fmt::format(" {} ", ref.func);
    // }
    return fmt::format("{:016b}{} | {} | c{}, d{}, l{} |", func, (xorable?'x':' '), F2STR[(US)last_func], cost, depth, lvl);
    // return fmt::format("{:04x}{} | {} | {}, {} | ", func, (xorable?'x':' '), F2STR[(US)last_func], cost, depth);
    // return "Func: " + std::to_string(func) + "|Last: " + std::to_string(last_func) + "|Cost: " + std::to_string(cost) + "|Depth: " + std::to_string(depth) + "|X: " + std::to_string(xorable);
  }
  
  std::string genlib_eqn(phmap::flat_hash_map<ULL, Node> & nodemap, std::vector<UI> & pis) const
  {
    // if (DEBUG) {fmt::print("\t\tAccessing genlib_eqn: {}\n", to_str());}
    if (last_func == fPI)
    {
      if (std::find(pis.begin(), pis.end(), func) == pis.end())
      {
        pis.push_back(func);
      }
      return PI2LETTER[func];
    }
    else if (last_func == fDFF)
    {
      assert(parent_hashes.size() == 1);
      return fmt::format("{}", nodemap[parent_hashes.back()].genlib_eqn(nodemap, pis));
    }
    else if (last_func == fNOT)
    {
      assert(parent_hashes.size() == 1);
      return fmt::format("!{}", nodemap[parent_hashes.back()].genlib_eqn(nodemap, pis));
    }
    else if (last_func == fCB)
    {
      assert(parent_hashes.size() == 2);
      return fmt::format("({0}+{1})", nodemap[parent_hashes.front()].genlib_eqn(nodemap, pis), nodemap[parent_hashes.back()].genlib_eqn(nodemap, pis));
    }
    else if (last_func == fOR || last_func == fMERGE)
    {
      assert(parent_hashes.size() == 2);
      return fmt::format("({0}|{1})", nodemap[parent_hashes.front()].genlib_eqn(nodemap, pis), nodemap[parent_hashes.back()].genlib_eqn(nodemap, pis));
    }
    else if (last_func == fAND)
    {
      assert(parent_hashes.size() == 2);
      return fmt::format("({0}&{1})", nodemap[parent_hashes.front()].genlib_eqn(nodemap, pis), nodemap[parent_hashes.back()].genlib_eqn(nodemap, pis));
    }
    else if (last_func == fXOR)
    {
      assert(parent_hashes.size() == 2);
      // return fmt::format("(!({0})*({1})+({0})*!({1}))", nodemap[parent_hashes.front()].genlib_eqn(nodemap, pis), nodemap[parent_hashes.back()].genlib_eqn(nodemap, pis));
      return fmt::format("({0}^{1})", nodemap[parent_hashes.front()].genlib_eqn(nodemap, pis), nodemap[parent_hashes.back()].genlib_eqn(nodemap, pis));
    }
    else
    {
      if (DEBUG) {fmt::print("Unsupported function {}", to_str());}
      return "";
    }
  }

  std::string to_genlib(phmap::flat_hash_map<ULL, Node> & nodemap, const std::vector<UI> & levels, const std::vector<UI> PI_funcs) const
  {
    std::vector<UI> pis;
    std::string str = fmt::format("GATE 0x{:04x}_{} {} O={};\n", func, fmt::join(levels, ""), cost,  genlib_eqn(nodemap, pis));
    for (auto & pi : pis)
    {
      UL idx = std::find(PI_funcs.begin(), PI_funcs.end(), pi) - PI_funcs.begin();
      auto true_lvl = (depth + 1) / 3;
      std::string line = fmt::format("\tPIN {} {} {} {} {:d} {:0.3f} {:d} {:0.3f}\n", 
      PI2LETTER[pi], GENLIB_PHASE, GENLIB_INPUT_LOAD, GENLIB_MAX_LOAD, 
      true_lvl - levels[idx], GENLIB_RISE_FANOUT_DELAY, true_lvl - levels[idx], GENLIB_FALL_FANOUT_DELAY);
      str.append(line);
      // if (DEBUG) {fmt::print(line);}
    }
    // if (DEBUG) {fmt::print(str);}
    return str;
  }

  // std::string to_genlib(phmap::flat_hash_map<ULL, Node> & nodemap, const std::vector<UI> & levels, phmap::flat_hash_map<UI, std::string> pi2symbol) const
  // {
  //   std::vector<UI> pis;
  //   std::string str = fmt::format("GATE 0x{:04x}_{} {} O={};\n", func, fmt::join(levels, ""), cost,  genlib_eqn(nodemap, pis));
  //   for (auto & pi : pis)
  //   {
  //     UL idx = std::find(PI_funcs.begin(), PI_funcs.end(), pi) - PI_funcs.begin();
  //     auto true_lvl = (depth + 1) / 3;
  //     std::string line = fmt::format("\tPIN {} {} {} {} {:d} {:0.3f} {:d} {:0.3f}\n", 
  //     PI2LETTER[pi], GENLIB_PHASE, GENLIB_INPUT_LOAD, GENLIB_MAX_LOAD, 
  //     true_lvl - levels[idx], GENLIB_RISE_FANOUT_DELAY, true_lvl - levels[idx], GENLIB_FALL_FANOUT_DELAY);
  //     str.append(line);
  //     // if (DEBUG) {fmt::print(line);}
  //   }
  //   // if (DEBUG) {fmt::print(str);}
  //   return str;
  // }

  UI node_cost(phmap::flat_hash_map<ULL, Node> & nodemap, const std::array<UI,12> cost_map) const
  {
    std::vector<ULL> stack { hash };
    UI gate_cost = 0;
    // fmt::print("Hash  inside {}\n", hash);
    // fmt::print("Node  inside {}\n", to_str());

    phmap::flat_hash_map<ULL, UI> ct_spl;
    std::unordered_set<ULL> non_splittable_nodes;

    // fmt::print("\t\t\tInitial cost {}\n", gate_cost);
    while (!stack.empty())
    {
      ULL n_hash = stack.back();
      stack.pop_back();
      Node& n = nodemap[n_hash]; //.at(n_hash); // [n_hash];
      // fmt::print("\t\t\tProcessing node {}\n", n.to_str());
      ct_spl[n_hash]++;
      if (ct_spl[n_hash] == 1)
      {
        gate_cost += cost_map[n.last_func];
        // fmt::print("\t\t\tFirst time, adding the cost of {} ({}). Total cost is {}\n", F2STR[n.last_func], cost_map[n.last_func], gate_cost);
        for (const auto& p_hash : n.parent_hashes)
        {
          stack.push_back(p_hash); // Node& p = nodemap[p_hash]; fmt::print("\t\t\tAdding parent node {} to stack\n", p.to_str());
        }
        if (n.last_func == fAND || n.last_func == fOR)
        {
          for (const auto& p_hash : n.parent_hashes)
          {
            non_splittable_nodes.emplace(p_hash);
          }
        }
      }
      else
      {
        gate_cost += cost_map[fSPL]; // 
        // fmt::print("\t\t\tNot first time, adding the cost of SPL ({}). Total cost is {}\n", cost_map[n.last_func], gate_cost);
      }
    }

    for (const ULL & n_hash : non_splittable_nodes)
    {
      auto count = ct_spl[n_hash];
      if (ct_spl[n_hash] > 1)
      {
        Node& n = nodemap[n_hash]; //(n_hash); //[n_hash];
        auto last_func = n.last_func;
        gate_cost += cost_map[last_func] * (count - 1); // duplicate a gate
        if (last_func == fXOR)
        {
          gate_cost += (count - 1) * cost_map[fSPL]; // need to do twice more splittings for duplicating an XOR gate
        }
      }
    }
    return gate_cost;
  }

  int recalculate_cost(phmap::flat_hash_map<ULL, Node> & nodemap, const std::array<UI,12> cost_map, std::vector<ULL> seen = {}) const
  {
    if (std::find(seen.begin(), seen.end(), hash) != seen.end())
    {
      return 0;
    }
    else
    {
      seen.push_back(hash);
    }
    if (last_func == fPI)
    {
      return 0;
    }
    else if (last_func == fDFF || last_func == fNOT)
    {
      return cost_map[last_func] + nodemap[parent_hashes.back()].recalculate_cost(nodemap, cost_map, seen);
    }
    else if (last_func == fCB || last_func == fOR || last_func == fAND || last_func == fMERGE || last_func == fXOR)
    {
      return cost_map[last_func] + nodemap[parent_hashes.front()].recalculate_cost(nodemap, cost_map, seen) + nodemap[parent_hashes.back()].recalculate_cost(nodemap, cost_map, seen);
    }
    else
    {
      return 9999;
    }
  }

  int number_of_dff(phmap::flat_hash_map<ULL, Node> & nodemap) const
  {
    // DFF is not redundant if 
    //    1. it is followed by SA gate
    //    2. it is followed by an XOR gate and preceded by a non-xorable node 
    
    std::vector<ULL> stack { hash };

    phmap::flat_hash_map<ULL, std::vector<ULL>> fanouts;
    phmap::flat_hash_map<ULL, std::vector<ULL>> fanins;

    while (!stack.empty())
    {
      ULL nhash = stack.back();
      stack.pop_back();
      Node & n = nodemap[nhash];

      for (auto phash : n.parent_hashes)
      {
        stack.push_back( phash );
        fanouts[phash].push_back(nhash);
        fanins[nhash].push_back(phash);
      }
    }

    int nPB_DFF = 0;
    for (auto const & [nhash, fo] : fanouts)
    {
      Node & n = nodemap[nhash];
      if (n.last_func != fDFF)
      {
        continue;
      }
      ULL pred_hash = fanins[nhash].back();
      Node & pred = nodemap[pred_hash];

      // check if followed by XOR
      if (!pred.xorable)
      {
        bool followed_by_xor = false;
        for (auto succ_hash : fo)
        {
          Node & succ = nodemap[succ_hash];
          if (succ.last_func == fXOR)
          {
            followed_by_xor = true;
            break;
          }
        }
        if (!followed_by_xor)
        {
          nPB_DFF++;
        }
      }
      else //DFF is not xorizing, check if there are any non-SA successors
      {
        bool followed_by_non_SA = true;
        for (auto succ_hash : fo)
        {
          Node & succ = nodemap[succ_hash];
          if (succ.last_func != fAND && succ.last_func != fOR)
          {
            followed_by_non_SA = false;
            break;
          }
        }
        if (!followed_by_non_SA)
        {
          nPB_DFF++;
        }
      }
    }

    return nPB_DFF;
  }

  std::vector<UI> structural_hash(phmap::flat_hash_map<ULL, Node> & nodemap) const
  {
    std::vector<ULL> stack { hash };
    std::vector<UI> nums;

    while (!stack.empty())
    {
      ULL nhash = stack.back();
      stack.pop_back();
      Node & n = nodemap[nhash];

      nums.push_back(n.last_func);
      nums.push_back(n.cost);
      
      std::vector<uint64_t> phashes = n.parent_hashes;
      std::sort(phashes.begin(), phashes.end(), 
        [&](const uint64_t& a, const uint64_t& b) 
        {
          if (nodemap[hash].last_func < nodemap[hash].last_func)
          {
            return true;
          }
          else if (nodemap[hash].last_func > nodemap[hash].last_func)
          {
            return false;
          }
          else if (nodemap[hash].last_func == nodemap[hash].last_func)
          {
            return (nodemap[hash].cost < nodemap[hash].cost);
          }

        }
      );
      stack.insert(stack.end(), phashes.begin(), phashes.end());
    }
    return nums;
  }

  void topo_sort(phmap::flat_hash_map<ULL, Node> & nodemap, std::vector<ULL>& order, bool verbose = false) const
  {
    // fmt::print("\t\tTraversing {0}\n", 
    //   hash
    // );
    for (const ULL phash : parent_hashes)
    {
      if (std::find(order.begin(), order.end(), phash) != order.end())
      {
        continue;
      }
      Node & p = nodemap[phash];
      p.topo_sort(nodemap, order, verbose);
    }
    order.push_back(hash);
    // fmt::print("Node {0}, func {1}, last_func {2}, parents: \n\t{3}", 
    //   hash,
    //   func,
    //   F2STR[last_func],
    //   fmt::join(parent_hashes, "\n\t")
    // );
  }

  std::string to_stack(phmap::flat_hash_map<ULL, Node> & nodemap, const phmap::flat_hash_map<uint16_t, uint16_t> & pi_map = dummy_map) const
  {
    if (last_func == fPI)
    {
      return PI2LETTER[func];
    }
    else if (last_func == fDFF)
    {
      assert(parent_hashes.size() == 1);
      return fmt::format("DFF({})", nodemap[parent_hashes.back()].to_stack(nodemap, pi_map));
    }
    else if (last_func == fNOT)
    {
      assert(parent_hashes.size() == 1);
      return fmt::format("NOT({})", nodemap[parent_hashes.back()].to_stack(nodemap, pi_map));
    }
    else if (last_func == fCB || last_func == fMERGE)
    {
      assert(parent_hashes.size() == 2);
      return fmt::format("CB({0}, {1})", nodemap[parent_hashes.front()].to_stack(nodemap, pi_map), nodemap[parent_hashes.back()].to_stack(nodemap, pi_map));
    }
    else if ( last_func == fOR )
    {
      assert(parent_hashes.size() == 2);
      return fmt::format("OR({0}, {1})", nodemap[parent_hashes.front()].to_stack(nodemap, pi_map), nodemap[parent_hashes.back()].to_stack(nodemap, pi_map));
    }
    else if (last_func == fAND)
    {
      assert(parent_hashes.size() == 2);
      return fmt::format("AND({0}, {1})", nodemap[parent_hashes.front()].to_stack(nodemap, pi_map), nodemap[parent_hashes.back()].to_stack(nodemap, pi_map));
    }
    else if (last_func == fXOR)
    {
      assert(parent_hashes.size() == 2);
      return fmt::format("XOR({0}, {1})", nodemap[parent_hashes.front()].to_stack(nodemap, pi_map), nodemap[parent_hashes.back()].to_stack(nodemap, pi_map));
    }
    else
    {
      if (DEBUG) {fmt::print("Unsupported function {}", to_str());}
      return "";
    }
  }

  auto process_nodes(phmap::flat_hash_map<ULL, Node> & GNM)
  {
    std::vector<std::pair<ULL,US>> stack { std::make_pair(hash, 0) };
    std::vector<ULL> seen;
    std::array<US,NUM_VARS> delays {{ 0xFF, 0xFF, 0xFF, 0xFF }};
    bool valid = true;
    while (!stack.empty())
    {
      auto [nhash, nlvl] = stack.back();
      stack.pop_back();
      
      if (VECTOR_CONTAINS(seen, nhash))
      {
        continue;
      }
      else
      {
        seen.push_back( nhash );
      }

      Node& n = GNM[nhash];
      if (n.last_func == fNOFUNC) 
      {
        valid = false;
        continue;
      }
      else if (n.last_func == fPI)
      {
        switch (n.func)
        {
          case 0x5555: delays[0] = nlvl; break;
          case 0x3333: delays[1] = nlvl; break;
          case 0x0F0F: delays[2] = nlvl; break;
          case 0x00FF: delays[3] = nlvl; break;
        }
        continue;
      }

      US new_lvl = nlvl + ((n.last_func == fDFF) | (n.last_func == fNOT) | (n.last_func == fXOR));

      for (ULL phash : n.parent_hashes)
      {
        stack.push_back( std::make_pair(phash, new_lvl) );
      }
    }

  }

  std::tuple<bool, UI> redundancy_check(phmap::flat_hash_map<ULL, Node> & GNM)
  {
    // std::vector<UI>   pi_funcs   {0x00FF, 0x0F0F, 0x3333, 0x5555};
    std::array<UI,  4> pi_funcs   {0x5555, 0x3333, 0x0F0F, 0x00FF};
    std::array<bool,4> has_dff    {false, false, false, false};
    std::array<bool,4> has_other  {false, false, false, false};
    std::array<bool,4> is_reached   {false, false, false, false};
    
    std::vector<ULL> stack = parent_hashes;
    if (DEBUG) {fmt::print("\tAnalyzing stack for node {}\n", to_str());}
    while (!stack.empty())
    {
      ULL n_hash = stack.back();
      stack.pop_back();
      Node& n = GNM[n_hash];
      if (n.last_func == fNOFUNC) continue;
      if (DEBUG) {fmt::print("\t\tAnalyzing node: {}\n", n.to_str());}
      auto it = std::find(pi_funcs.begin(), pi_funcs.end(), n.func);
      if (it != pi_funcs.end()) // if the function is a PI
      {
        auto idx = it - pi_funcs.begin();
        if (DEBUG) {fmt::print("\t\tFound PI at idx {} for {:04x}\n", idx, n.func);}
        // if (DEBUG) {fmt::print("\t\tLast_func == {}\n", n.last_func);}
        // if (DEBUG) {fmt::print("\t\fDFF == {}\n", fDFF);}
        if (n.last_func == fDFF) 
        {
          has_dff[idx] = true;
        }
        else if (n.last_func == fPI)
        {
          is_reached[idx] = true;
        }
        else 
        {
          has_other[idx] = true;
        }
        if (DEBUG) {fmt::print("\t\t\tNew has_dff:\t{}\n", fmt::join(has_dff, "\t"));}
        if (DEBUG) {fmt::print("\t\t\tNew has_other:\t{}\n", fmt::join(has_other, "\t"));}
        if (DEBUG) {fmt::print("\t\t\tNew is_reached:\t{}\n", fmt::join(is_reached, "\t"));}
      }
      /*
        if (n.parent_hashes.size() > 1)
        {
          bool all_dff = true;
          for (auto phash : n.parent_hashes)
          {
            Node & p = GNM.at(phash);
            if (p.last_func != fDFF) 
            {
              all_dff = false;
              break;
            }
          }
          if (all_dff)
          {
            return std::make_tuple(false, 0);
          }
        }
      */
      stack.insert(stack.end(), n.parent_hashes.begin(), n.parent_hashes.end());
    }   

    UI support_size = NUM_VARS;

    if (DEBUG) {fmt::print("\t\tAnalyzing vectors\n");}
    for (auto i = 0u; i < NUM_VARS; ++i)
    {
      if (DEBUG) {fmt::print("\t\tPI : {}\t func:{:04x} \t has_dff: {} | has_other: {}| is_reached: {}\n", i, pi_funcs[i], has_dff[i],  has_other[i], is_reached[i]);   }
      if (has_other[i])
      {
        if (DEBUG) {fmt::print("\t\t\t OK PI\n");}
        continue;
      }
      else if ( has_dff[i] ) 
      {
        assert(is_reached[i]); 
        if (DEBUG) {fmt::print("\t\t\t Violating PI\n");}
        return std::make_tuple(false, 0);
      }
      else if ( ~has_dff[i] ) 
      {       
        if (is_reached[i]) 
        {
          if (DEBUG) {fmt::print("\t\t\t OK PI\n");}
          continue;
        }
        else
        {
          if (DEBUG) {fmt::print("\t\t\t Redundant PI\n");}
          support_size--; 
          continue;
        }
      }
    }
    return std::make_tuple(true, support_size);
  }

private:

};

namespace std {
  template<> struct hash<Node> {
    std::size_t operator()(const Node& node) const {
      return std::hash<ULL>{}(node.hash);
    }
  };
}


std::vector<UI> gen_pi_func(UI nvars)
{
  if (nvars == 1)
  {
    return {1};
  }
  UI power = (1 << (1 << (nvars - 1)));
  UI factor = power + 1;
  std::vector<UI> out = {power - 1};
  for (auto n : gen_pi_func(nvars - 1))
  {
    out.push_back(factor * n);
  }
  return out;
}

template <std::size_t N>
void apply_permutation(const std::vector<int>& perm, std::array<int, N>& delay) 
{
  std::array<int, N> temp;
  for (std::size_t i = 0; i < N; ++i) 
  {
    temp[i] = delay[perm[i]];
  }
  delay = temp;
}

/* Returns the functions of PIs and their correponding delays*/
std::tuple<
  std::array<uint8_t, NUM_VARS>, 
  std::vector<uint8_t>, 
  bool
  > get_delays(
    uint64_t hash, 
    phmap::flat_hash_map<ULL, Node> & hashmap
  )
{
  std::vector<std::pair<uint64_t, uint8_t>> stack;
  stack.emplace_back(hash, 0);
  std::vector<uint64_t> seen;
  std::array<uint8_t, NUM_VARS> delays;
  std::fill(delays.begin(), delays.end(), 0xFF);
  std::vector<uint8_t> support;

  Node & root = hashmap[hash];
  /* Special case for constants */
  if (root.func == 0xFFFF || root.func == 0x0000)
  {
    return std::make_tuple( delays, support, false );
  }

  while (!stack.empty())
  {
    auto [h, delay] = stack.back();
    Node & n = hashmap[h];
    // fmt::print("\tStack: analyzing {}: {}\n", h, n.to_str());
    stack.pop_back();
    if (std::find( seen.begin(), seen.end(), h) != seen.end() )
    {
      continue;
    }
    seen.push_back(h);

    // Node & n = hashmap[h];
    /* Invalid node encountered, discard the TT */
    if (n.last_func == fNOFUNC)
    {
      return std::make_tuple( delays, support, false );
    }
    /* Found PI, update delay and support */
    else if (n.last_func == fPI)
    {
      uint8_t idx = PIFUNC2IDX.at(n.func);
      delays[idx] = delay;
      support.push_back(idx);
      continue;
    }

    /* Regular node, add parents to stack. Increment delay if the last element is clocked */
    uint8_t new_delay = delay + (n.last_func == fDFF || n.last_func == fNOT || n.last_func == fXOR);
    for (uint64_t phash : n.parent_hashes)    
    {
      stack.emplace_back( phash, new_delay );
    }
  }
  /* Success, return the result */
  return std::make_tuple( delays, support, true );
}

std::tuple<TT4, 
  std::vector<uint8_t>, 
  uint8_t,
  std::string,
  bool> process_node(
    const uint64_t hash, 
    phmap::flat_hash_map<ULL, Node> & hashmap)
{
  auto [pi_delays, support, status] = get_delays(hash, hashmap);
  if (!status) return std::make_tuple( TT4{}, std::vector<uint8_t>{}, 0, "", false);

  Node & n = hashmap[hash];
  TT4 tt;
  tt._bits = n.func;
  const auto num_vars = tt.num_vars();

  /* Mask for cases when support size < 4 */
  TT4 mask;
  mask._bits = (1 << (1 << support.size())) - 1;

  auto t1 = tt;
  auto tmin = t1;

  phmap::flat_hash_map<uint16_t, uint8_t> pi_map;
  for (auto idx : support)
  {
    uint16_t func = PI_WORDS[idx];
    pi_map.emplace(PI_WORDS[idx], idx);
  }

  // std::vector<uint8_t> pi_order(support.size());
  // std::iota(pi_order.begin(), pi_order.end(), 0);

  std::string chars = "abcd";
  std::string best_chars = chars;
  auto new_pi_delays = pi_delays; // to store the swapped delays
  // auto new_pi_order = pi_order; // to store the swapped delays
  auto new_pi_map = pi_map;
  for (auto j = support.size(); j < NUM_VARS; ++j)
  {
    if ( kitty::has_var( t1, j) )
    {
      tmin._bits = 0xFFFF;
      break;
    }
  }
  const auto& swaps = kitty::detail::swaps[num_vars - 2u];
  int best_swap = -1;
  for ( std::size_t i = 0; i < swaps.size(); ++i )
  {
    const auto pos = swaps[i];
    kitty::swap_adjacent_inplace( t1, pos );
    std::swap(chars[pos], chars[pos+1]);
    std::swap(pi_delays[pos], pi_delays[pos+1]);
    // std::swap(pi_order[pos], pi_order[pos+1]);


    /* Check if TT is min_base */
    bool acceptable = true;
    for (auto j = support.size(); j < NUM_VARS; ++j)
    {
      if ( kitty::has_var( t1, j) )
      {
        acceptable = false;
        break;
      }
    }
    
    /* If TT is lexicographically smallest and is min_base */
    if ( ( (t1 & mask) < (tmin & mask) ) && acceptable)
    {
      best_swap = static_cast<int>( i );
      tmin = t1;
      new_pi_delays = pi_delays;
      best_chars = chars;
      // new_pi_order = pi_order;
    }
  }
  std::vector<uint8_t> new_pi_delays_vector;
  new_pi_delays_vector.insert(
    new_pi_delays_vector.end(), 
    new_pi_delays.begin(), 
    new_pi_delays.begin() + support.size()
  );
  return std::make_tuple( tmin, new_pi_delays_vector, support.size(), best_chars, true );
}


#pragma region write_output

void write_csv_gnm(const phmap::flat_hash_map<ULL, Node>& gnm, const std::string& filename) {
  // Open output file
  std::ofstream outfile(filename);

  // Write header row to CSV file
  outfile << "Hash,Func,Last Func,Cost,Depth,Xorable,Parent Hashes,Lvl" << std::endl;

  // Write data to CSV file
  for (const auto& [hash, n] : gnm) {
    std::string str = fmt::format("{0},{1},{2},{3},{4},{5:d},{6},{7}", 
                    hash, n.func, n.last_func, n.cost, n.depth, n.xorable, 
                    fmt::join(n.parent_hashes, "|"), n.lvl);
    outfile << str << std::endl;
  }

  // Close output file
  outfile.close();
}

void write_csv_arr(const std::array<ULL, NUM_TT>& arr_hashes, const std::string& filename) {
  // Open output file
  std::ofstream outfile(filename);

  // Write header row to CSV file
  outfile << "Hash" << std::endl;

  // Write data to CSV file
  for (const auto& hash : arr_hashes) {
    outfile << fmt::format("{0}", hash) << std::endl;
  }
  // Close output file
  outfile.close();
}

phmap::flat_hash_map<ULL, Node> read_csv_gnm(const std::string& filename, const uint32_t cost_cutoff = 512u, const uint32_t depth_cutoff = 16u) 
{
  // Open input file
  std::ifstream infile(filename);
  
  phmap::flat_hash_map<ULL, Node> gnm;

  // Parse CSV file and populate GNM variable
  std::string line;
  std::getline(infile, line);  // skip header row
  while (std::getline(infile, line)) 
  {
    std::stringstream ss;
    ss.str(line);
    std::string field;
    Node node;
    std::getline(ss, field, ',');
    std::stringstream(field) >> node.hash;
    std::getline(ss, field, ',');
    std::stringstream(field) >> node.func;
    std::getline(ss, field, ',');
    std::stringstream(field) >> node.last_func;
    if (node.last_func > 20) { continue; }
    std::getline(ss, field, ',');
    std::stringstream(field) >> node.cost;
    if (node.cost > cost_cutoff) { continue; }
    std::getline(ss, field, ',');
    std::stringstream(field) >> node.depth;
    if (node.depth > depth_cutoff) { continue; }
    std::getline(ss, field, ',');
    std::stringstream(field) >> node.xorable;
    std::getline(ss, field, ',');
    std::stringstream parent_hashes_ss(field);
    while (std::getline(parent_hashes_ss, field, '|')) 
    {
      ULL parent_hash;
      std::stringstream(field) >> parent_hash;
      node.parent_hashes.push_back(parent_hash);
    }
    std::getline(ss, field, ',');
    std::stringstream(field) >> node.lvl;
    gnm.emplace(node.hash, node);
  }

  // Close input file
  infile.close();
  return gnm;
}

std::array<ULL, NUM_TT> read_csv_arr(const std::string& filename) {
  // Open input file
  std::ifstream infile(filename);
  std::array<ULL, NUM_TT> arr_hashes;

  // Parse CSV file and populate array
  std::string line;
  std::getline(infile, line);  // skip header row
  UI index = 0;
  while (std::getline(infile, line) && index < NUM_TT) {
    std::stringstream ss(line);
    std::string field;
    std::getline(ss, field, ',');
    std::stringstream(field) >> arr_hashes[index++];
  }

  // Close input file
  infile.close();
  return arr_hashes;
}


bool is_good(ULL hash, Node & node, phmap::flat_hash_map<ULL, bool> & status, phmap::flat_hash_map<ULL, Node>& all_hashes)
{
  bool good = true;
  for (auto & p_hash : node.parent_hashes)
  {
    Node & p_node = all_hashes[p_hash];
    if (status.find(p_hash) == status.end()) // if a parent status is unknown
    {
      status[p_hash] = is_good(p_hash, p_node, status, all_hashes);
    }
    good &= status[p_hash];
  }
  return good;
}

phmap::flat_hash_map<ULL, bool> subset_of_pi(std::vector<ULL>& pi, phmap::flat_hash_map<ULL, Node>& all_hashes)
{
  phmap::flat_hash_map<ULL, bool> status;
  for (ULL hash : pi)
  {
    status[hash] = true;
  }

  for (auto & [hash, node] : all_hashes)
  {
    status[hash] = is_good(hash, node, status, all_hashes);
  }
  return status;
}

#pragma endregion

std::pair<std::string, bool> get_formula_old( TT4 & tt, int var_idx = -1, bool reverse = false, const std::vector<std::string> & varnames = {"a","b","c","d"})
{
  if (reverse)
  {
    // TT4 tt_original = tt;
    fmt::print("\tTT before: {}\n", kitty::to_hex(tt));
    uint16_t num_orig = tt._bits;
    uint16_t num_new = 0u;
    for (auto i = 0u; i < 16; ++i)
    {
      // kitty::copy_bit(tt_original, i, tt, 16 - i - 1 );
      num_new |= ((num_orig >> i) & 1) << (16 - i - 1);
    }
    tt._bits = num_new;
    fmt::print("\tTT after: {}\n", kitty::to_hex(tt));
  }
  if (kitty::is_const0(tt))
  {
    return std::make_pair("", false);
  }
  else if (kitty::is_const0(~tt))
  {
    return std::make_pair("", true);
  }
  else
  {
    var_idx++;
    TT4 tt0 = kitty::cofactor0(tt, var_idx);
    TT4 tt1 = kitty::cofactor1(tt, var_idx);

    auto [s0, valid0] = get_formula_old(tt0, var_idx, false, varnames);
    auto [s1, valid1] = get_formula_old(tt1, var_idx, false, varnames);

    std::string out0;
    std::string out1;
    if (valid0)
    {
      if (s0.size() > 0)
      {
        out0 = fmt::format("(!{}&{})", varnames[var_idx], s0);
      }
      else
      {
        out0 = fmt::format("(!{})", varnames[var_idx]);
      }
    }
    if (valid1)
    {
      if (s1.size() > 0)
      {
        out1 = fmt::format("({}&{})", varnames[var_idx], s1);
      }
      else
      {
        out1 = fmt::format("({})", varnames[var_idx]);
      }
    }

    if (out0.size() > 0)
    {
      if (out1.size() > 0)
      {
        return std::make_pair(fmt::format("({0}|{1})", out0, out1), true);
      }
      else
      {
        return std::make_pair(out0, true);
      }
    }
    else
    {
      if (out1.size() > 0)
      {
        return std::make_pair(out1, true);
      }
      else
      {
        return std::make_pair("", false);
      }
    }
  }
}

template <typename TT> 
TT reverse_TT(TT & tt)
{

  // fmt::print("\tTT before: {}\n", kitty::to_hex(tt));
  uint16_t num_orig = tt._bits;
  uint16_t num_new = 0u;
  int num_vars = tt.num_vars();
  for (auto i = 0u; i < (1 << (1 << num_vars)); ++i)
  {
    num_new |= ((num_orig >> i) & 1) << ((1 << num_vars) - i - 1);
  }
  tt._bits = num_new;
  // fmt::print("\tTT after: {}\n", kitty::to_hex(tt));
  return tt;
}

std::string get_formula( TT4 _tt, const bool reverse = false, const std::vector<std::string> & varnames = {"a","b","c","d"})
{
  std::vector<std::pair<TT4, std::vector<std::string>>> stack;
  if (reverse)
  { 
    stack.emplace_back(reverse_TT(_tt), std::vector<std::string>{});
  }
  else
  { 
    stack.emplace_back(       _tt , std::vector<std::string>{});
  };

  std::vector<std::vector<std::string>> minterms;

  while (!stack.empty())
  {
    auto [tt, minterm] = stack.back();
    stack.pop_back();

    if (kitty::is_const0(tt))
    {
      continue;
    }
    else if (kitty::is_const0(~tt))
    {
      
      minterms.push_back(minterm);
      continue;
    }
    else
    {
      int var_idx = minterm.size();
      TT4 tt0 = kitty::cofactor0(tt, var_idx);
      std::vector<std::string> minterm0 = minterm;
      minterm0.push_back("!" + varnames[var_idx]);
      // minterm0.push_back(fmt::format("(0xFFFF^{})",varnames[var_idx]));
      stack.emplace_back( tt0, minterm0 );
      // fmt::print("Adding minterm : {}\n", fmt::join(minterm0,", "));

      TT4 tt1 = kitty::cofactor1(tt, var_idx);
      std::vector<std::string> minterm1 = minterm;
      minterm1.push_back(    varnames[var_idx]);
      stack.emplace_back( tt1, minterm1 );
      // fmt::print("Adding minterm : {}\n", fmt::join(minterm1,", "));
    }
  }
  
  std::vector<std::string> out;
  for (std::vector<std::string> & minterm : minterms )
  {
    if (minterm.size() > 1)
    {
      out.push_back(fmt::format("({})", fmt::join(minterm, "&")));
    }
    else
    {
      out.push_back(fmt::format( "{}" , fmt::join(minterm, "&")));
    }
  }
  
  TT4 tt_test;
  std::string expression = fmt::format( "{}" , fmt::join(out, "|"));
  kitty::create_from_formula(tt_test, expression, varnames);
  assert(tt_test == _tt);
  return expression;
}

// phmap::flat_hash_map<ULL, Node> mergeMaps(const std::vector<phmap::flat_hash_map<ULL, Node>>& maps) 
// {
//     phmap::flat_hash_map<ULL, Node> mergedMap;

//     for (const auto& map : maps) {
//         mergedMap.insert(map.begin(), map.end());
//     }

//     return mergedMap;
// }

phmap::flat_hash_map<ULL, Node> mergeMaps(std::vector<phmap::flat_hash_map<ULL, Node>>& maps) 
{
  // phmap::flat_hash_map<ULL, Node> mergedMap;

  // for (const auto& map : maps) 
  for (auto it = maps.rbegin() - 1; it < maps.rend(); it++) 
  {
    auto & other_map = *it;
    maps[0].insert(other_map.begin(), other_map.end());
  }

  return maps[0];
}

phmap::flat_hash_map<ULL, Node> read_global_gnm(
  const std::vector<std::vector<UI>> sets_of_levels,
  const std::string file_prefix,
  const std::array<int, 12> costs = COSTS_CONNECT_CONSERVATIVE
)
{
  std::vector<phmap::flat_hash_map<ULL, Node>> GNM_vector(sets_of_levels.size());
  std::vector<std::thread> threads;

  // Launch threads to read the CSV files in parallel
  auto i = 0u;
  for (const std::vector<UI> & levels : sets_of_levels)
  {
    threads.emplace_back([i, &GNM_vector, &file_prefix, &levels, &costs]() 
    {
      GNM_vector[i] = read_csv_gnm(fmt::format("{}_{}_gnm.csv", file_prefix, fmt::join(levels, "")));
      fmt::print("Imported GNM {0} with {1} nodes\n", fmt::join(levels, ""), GNM_vector[i].size());
      // for (auto & [hash, node] : GNM_vector[i])
      // {
      //   // fmt::print("Imported node {}\n", node.to_str());
      //   UI node_cost = node.node_cost(GNM_vector[i], costs);

      //   node.cost = node_cost;
      //   // GNM_global[hash].cost = node_cost;
      //   // GNM_global[hash] = node;
      // }
    });
    i++;
  }
  // Wait for all threads to finish
  for (auto& thread : threads) 
  {
      thread.join();
  }

  phmap::flat_hash_map<ULL, Node> GNM_global = mergeMaps(GNM_vector);
  return GNM_global;
}

void SaveToFile(const phmap::flat_hash_map<uint64_t, Node>& map, const std::string& filepath_prefix, uint64_t max_size = 1024 * 1024 * 50) 
{
  fmt::print("Starting the dumping process\n");

  auto it = map.begin();
  int chunk_id = 0;

  while (it != map.end())
  {
    std::string filepath = fmt::format("{}_{}.dat", filepath_prefix, chunk_id);
    fmt::print("Opening {} for writing\n", filepath);
    std::ofstream outfile(filepath, std::ios::binary);
    if (!outfile) 
    {
      std::cerr << "Failed to open file for writing" << std::endl;
      return;
    }
    // if (chunk_id == 0)
    // {
    //   // Write out the size of the map first
    //   size_t map_size = map.size();
    //   outfile.write(reinterpret_cast<const char*>(&map_size), sizeof(map_size));
    //   fmt::print("Recorded map size = {}\n", map_size);
    // }
    
    while ( (outfile.tellp() < max_size) && (it != map.end()) )
    {
      auto & [hash, node] = *it;
      outfile.write(reinterpret_cast<const char*>(&node.hash), sizeof(ULL));

      outfile.write(reinterpret_cast<const char*>(&node.func), sizeof(UI));
      outfile.write(reinterpret_cast<const char*>(&node.last_func), sizeof(US));
      outfile.write(reinterpret_cast<const char*>(&node.cost), sizeof(UI));
      outfile.write(reinterpret_cast<const char*>(&node.depth), sizeof(UI));
      outfile.write(reinterpret_cast<const char*>(&node.xorable), sizeof(bool));
      
      // Write vector size and contents
      US vec_size = node.parent_hashes.size();
      outfile.write(reinterpret_cast<const char*>(&vec_size), sizeof(US));
      if (vec_size > 0) 
      {
        outfile.write(reinterpret_cast<const char*>(&node.parent_hashes[0]), vec_size * sizeof(ULL));
      }

      outfile.write(reinterpret_cast<const char*>(&node.lvl), sizeof(UI));
      it++;
    }
    // if the map is done, place a codeword to the end of the file
    // if (it == map.end())
    // {
    //   const ULL code = 0x1c0fb3767dc1ea0c;
    //   outfile.write(reinterpret_cast<const char*>(&code), sizeof(ULL));
    // }

    outfile.close();
    chunk_id++;
  } 
}

bool LoadFromFile(phmap::flat_hash_map<uint64_t, Node>& map, const std::string& filepath_prefix) 
{
  fmt::print("Starting the loading process\n");

  int chunk_id = 0;
  while (true)
  {
    std::string filepath = fmt::format("{}_{}.dat", filepath_prefix, chunk_id);
    std::ifstream infile(filepath, std::ios::binary);
    if (!infile) 
    {
      fmt::print("Reading complete\n", filepath);
      return true;
    }
    fmt::print("Reading from {}\n", filepath);
    
    // Determine the size of the file in bytes
    infile.seekg(0, std::ios::end);
    std::streampos fileSize = infile.tellg();
    infile.seekg(0, std::ios::beg);

    while (infile.tellg() < fileSize)
    {
      Node node;
      infile.read(reinterpret_cast<char*>(&node.hash), sizeof(ULL));
      // if (node.hash == 0x1c0fb3767dc1ea0c)
      // {
      //   // found the end of the map
      //   infile.close();
      //   return true;
      // }

      infile.read(reinterpret_cast<char*>(&node.func), sizeof(UI));
      infile.read(reinterpret_cast<char*>(&node.last_func), sizeof(US));
      infile.read(reinterpret_cast<char*>(&node.cost), sizeof(UI));
      infile.read(reinterpret_cast<char*>(&node.depth), sizeof(UI));
      infile.read(reinterpret_cast<char*>(&node.xorable), sizeof(bool));

      US vec_size = 0;
      infile.read(reinterpret_cast<char*>(&vec_size), sizeof(vec_size));
      node.parent_hashes.resize(vec_size);
      if (vec_size > 0) 
      {
        infile.read(reinterpret_cast<char*>(&node.parent_hashes[0]), vec_size * sizeof(ULL));
      }

      infile.read(reinterpret_cast<char*>(&node.lvl), sizeof(UI));
      // fmt::print("Loading: {}\n", node.to_str());
      map.emplace( std::move(node).hash, std::move(node) );
    }
    infile.close();
    chunk_id++;
  }
}

phmap::flat_hash_map<uint64_t, Node> LoadFromOneFile(const int chunk_id, const std::string & path) 
{
  phmap::flat_hash_map<uint64_t, Node> map;
  std::ifstream infile(path, std::ios::binary);
  fmt::print("Reading from {}\n", path);
  
  // Determine the size of the file in bytes
  infile.seekg(0, std::ios::end);
  std::streampos fileSize = infile.tellg();
  infile.seekg(0, std::ios::beg);

  while (infile.tellg() < fileSize)
  {
    Node node;
    infile.read(reinterpret_cast<char*>(&node.hash), sizeof(ULL));
    infile.read(reinterpret_cast<char*>(&node.func), sizeof(UI));
    infile.read(reinterpret_cast<char*>(&node.last_func), sizeof(US));
    infile.read(reinterpret_cast<char*>(&node.cost), sizeof(UI));
    infile.read(reinterpret_cast<char*>(&node.depth), sizeof(UI));
    infile.read(reinterpret_cast<char*>(&node.xorable), sizeof(bool));

    US vec_size = 0;
    infile.read(reinterpret_cast<char*>(&vec_size), sizeof(vec_size));
    node.parent_hashes.resize(vec_size);
    if (vec_size > 0) 
    {
      infile.read(reinterpret_cast<char*>(&node.parent_hashes[0]), vec_size * sizeof(ULL));
    }
    infile.read(reinterpret_cast<char*>(&node.lvl), sizeof(UI));
    map.emplace( std::move(node).hash, std::move(node) );
    fmt::print("[{}] Loading: {}\n", chunk_id,  node.to_str());
  }
  infile.close();
  return map;
}

void LoadVectorFromOneFile(std::vector<Node> & vec, const int chunk_id, const std::string & path) 
{
  std::ifstream infile(path, std::ios::binary);
  fmt::print("Reading from {}\n", path);
  
  // Determine the size of the file in bytes
  infile.seekg(0, std::ios::end);
  std::streampos fileSize = infile.tellg();
  infile.seekg(0, std::ios::beg);

  while (infile.tellg() < fileSize)
  {
    Node node;
    infile.read(reinterpret_cast<char*>(&node.hash), sizeof(ULL));
    infile.read(reinterpret_cast<char*>(&node.func), sizeof(UI));
    infile.read(reinterpret_cast<char*>(&node.last_func), sizeof(US));
    infile.read(reinterpret_cast<char*>(&node.cost), sizeof(UI));
    infile.read(reinterpret_cast<char*>(&node.depth), sizeof(UI));
    infile.read(reinterpret_cast<char*>(&node.xorable), sizeof(bool));

    US vec_size = 0;
    infile.read(reinterpret_cast<char*>(&vec_size), sizeof(vec_size));
    node.parent_hashes.resize(vec_size);
    if (vec_size > 0) 
    {
      infile.read(reinterpret_cast<char*>(&node.parent_hashes[0]), vec_size * sizeof(ULL));
    }
    infile.read(reinterpret_cast<char*>(&node.lvl), sizeof(UI));
    fmt::print("[{}] Loading: {}\n", chunk_id,  node.to_str());
    
    vec.push_back( node );
  }
  infile.close();
}

phmap::flat_hash_map<ULL, Node> ParallelLoadFromFiles(const std::string& filepath_prefix) 
{
  fmt::print("Starting the loading process\n");
  fmt::print("Building the list of files\n");

  int chunk_id = 0;
  do
  {
    std::string filepath = fmt::format("{}_{}.dat", filepath_prefix, chunk_id);
    std::ifstream infile(filepath, std::ios::binary);
    if (!infile) 
    {
      break;
    }
    chunk_id++;
  } while (true);

  std::vector<std::thread> threads;
  threads.reserve(chunk_id);
  // std::vector<phmap::flat_hash_map<ULL, Node>> maps( chunk_id );
  std::vector<std::vector<Node>> vecs;
  for (int i = 0u; i < chunk_id; ++i)
  {
    // fmt::print("Reading {}_{}.dat\n", filepath_prefix, i);
    threads.emplace_back([i, &vecs, &filepath_prefix]() 
    {
      std::string filepath = fmt::format("{}_{}.dat", filepath_prefix, i);
      LoadVectorFromOneFile( vecs[i], i, filepath );
      // fmt::print("[{}] found {} nodes\n", i, maps[i].size());
    });
  }
  for (auto& thread : threads) 
  {
      thread.join();
  }

  phmap::flat_hash_map<ULL, Node> global_map;
  for (std::vector<Node> & vec: vecs)
  {
    for (Node & node : vec)
    {
      global_map.emplace(node.hash, std::move(node));
    }
  }
  return global_map;
}


#if false
  // Save the map to a binary file
  bool SaveToFile(const std::string& filename, const std::vector<Node>& data) {
      std::ofstream file(filename, std::ios::binary | std::ios::out);
      if (!file.is_open()) {
          std::cerr << "Error: Unable to open file for writing: " << filename << std::endl;
          return false;
      }

      // Write the number of nodes
      size_t numNodes = data.size();
      file.write(reinterpret_cast<const char*>(&numNodes), sizeof(size_t));

      // Write the nodes
      for (const Node& node : data) {
          file.write(reinterpret_cast<const char*>(&node), sizeof(Node));
      }

      // Close the file
      file.close();
      return true;
  }

  // Load the map from a binary file
  bool LoadFromFile(const std::string& filename, std::vector<Node>& data) {
      std::ifstream file(filename, std::ios::binary | std::ios::in);
      if (!file.is_open()) {
          std::cerr << "Error: Unable to open file for reading: " << filename << std::endl;
          return false;
      }

      // Read the number of nodes
      size_t numNodes;
      file.read(reinterpret_cast<char*>(&numNodes), sizeof(size_t));

      // Resize the vector to hold the data
      data.resize(numNodes);

      // Read the nodes
      for (Node& node : data) {
          file.read(reinterpret_cast<char*>(&node), sizeof(Node));
      }

      // Close the file
      file.close();
      return true;
  }
#endif

#if false
  void DumpToFile(const phmap::flat_hash_map<uint64_t, Node>& map, const std::string& filename) {
    std::ofstream outfile(filename, std::ios::binary);

    if (!outfile) {
      std::cerr << "Error opening file for writing: " << filename << std::endl;
      return;
    }

    // Write the number of elements in the map
    uint64_t num_elements = map.size();
    outfile.write(reinterpret_cast<const char*>(&num_elements), sizeof(uint64_t));

    // Write each key-value pair to the file
    for (const auto& entry : map) {
      uint64_t key = entry.first;
      Node node = entry.second;

      outfile.write(reinterpret_cast<const char*>(&key), sizeof(uint64_t));
      outfile.write(reinterpret_cast<const char*>(&node), sizeof(Node));
    }

    outfile.close();
  }

    phmap::flat_hash_map<uint64_t, Node> LoadFromFile(const std::string& filename) {
      phmap::flat_hash_map<uint64_t, Node> map;
      std::ifstream infile(filename, std::ios::binary);

      if (!infile) {
        std::cerr << "Error opening file for reading: " << filename << std::endl;
        return map; // Return an empty map if there was an error
      }

      // Read the number of elements in the file
      uint64_t num_elements;
      infile.read(reinterpret_cast<char*>(&num_elements), sizeof(uint64_t));

      // Read each key-value pair from the file and insert it into the map
      for (uint64_t i = 0; i < num_elements; ++i) 
      {
        uint64_t key;
        Node node;

        infile.read(reinterpret_cast<char*>(&key), sizeof(uint64_t));
        infile.read(reinterpret_cast<char*>(&node), sizeof(Node));

        map[key] = node;
      }

      infile.close();
      return map;
    }
#endif