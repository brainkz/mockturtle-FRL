#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <fmt/format.h>
#include "mockturtle/algorithms/nodes.hpp"

constexpr uint8_t MAX_PER_FUNC = 5u;

struct Delay 
{
  std::vector<uint8_t> delays;
  Delay() = default;

  explicit Delay(const std::vector<uint8_t>& v) 
  {
    delays = v;
    // std::vector<uint8_t> delays;
    // std::copy(v.begin(), v.end(), delays.begin());
  }
  template <typename collection> 
  explicit Delay(const collection c) 
  {
    std::vector<uint8_t> delays;
    std::copy(c.begin(), c.end(), delays.begin());
  }

  bool strictly_better_than(const Delay& other) const 
  {
    if (delays.size() != other.delays.size())
    {
      return false;
    }
    bool strictly_smaller = false;
    for (size_t i = 0; i < delays.size(); ++i) 
    {
      if (delays[i] < other.delays[i]) 
      {
        strictly_smaller = true;
      } else if (delays[i] > other.delays[i]) 
      {
        return false;
      }
    }
    return strictly_smaller;
  }
  bool operator==(const Delay& other) const 
  {
    return delays == other.delays;
  }
  uint8_t size()
  {
    return delays.size();
  }
};

namespace std 
{
  template <>
  struct hash<Delay> 
  {
    size_t operator()(const Delay& d) const 
    {
        size_t h = 0;
        for (size_t i = 0; i < d.delays.size(); ++i) 
        {
        h = h * 31 + std::hash<uint8_t>()(d.delays[i]);
        }
        return h;
    }
  };
}

struct LibEntry
{
  std::string name;
  Delay delays;
  uint64_t hash;
  std::vector<uint8_t> pi_order;  
  /*
    PI order is [1,0,3] means that inputs [a,b,c] to a node in kLUT ntk correspond to PIs [0x3333, 0x5555, 0x00FF], respectively. 
  */
  std::vector<UI> struct_hash;

  // Copy constructor
  LibEntry(const LibEntry& other) :
    name(other.name),
    delays(other.delays),
    hash(other.hash),
    pi_order(other.pi_order), 
    struct_hash(other.struct_hash) {}

  // Explicit constructor
  explicit LibEntry
  (
    const std::string _name,
    const Delay& _delays, 
    const uint64_t _hash, 
    const std::vector<uint8_t>& _pi_order,
    std::unordered_map<uint64_t, Node> & nodemap
  ) :
    name(_name),
    delays(_delays),
    hash(_hash),
    pi_order(_pi_order) {
      struct_hash = nodemap[hash].structural_hash(nodemap);
    }

  explicit LibEntry(const std::string& entryStr) 
  {
    std::istringstream iss(entryStr);
    std::string token;

    std::getline(iss, name, ',');
    while (std::getline(iss, token, '|')) 
    {
        delays.delays.push_back(static_cast<uint8_t>(std::stoi(token)));
    }
    std::getline(iss, token, ',');
    hash = std::stoull(token);
    while (std::getline(iss, token, '|')) 
    {
        pi_order.push_back(static_cast<uint8_t>(std::stoi(token)));
    }
  }

  // Assignment operator
  LibEntry& operator=(const LibEntry& other) 
  {
    name = other.name;
    delays = other.delays;
    hash = other.hash;
    pi_order = other.pi_order;
    struct_hash = other.struct_hash;
    return *this;
  };
  UI cost(std::unordered_map<uint64_t, Node> & nodemap) const
  {
    return nodemap[hash].cost;
  }
  // UI cost(std::unordered_map<uint64_t, Node> & nodemap) const
  // {
  //   return nodemap[hash].cost;
  // }

  std::string write() const
  {
    return fmt::format(
      "{0},{1},{2},{3}\n",
      name,
      fmt::join(delays.delays, "|"),
      hash,
      fmt::join(pi_order, "|")
    );
  }
};

void writeLibEntryToCSV(const std::vector<LibEntry>& entries, const std::string& filename) 
{
  std::ofstream outFile(filename);
  if (!outFile) 
  {
    std::cerr << "Error opening file for writing: " << filename << std::endl;
    return;
  }

  for (const LibEntry& entry : entries) 
  {
    outFile << entry.write();
  }
  outFile.close();
}

std::vector<LibEntry> readLibEntryFromCSV(const std::string& filename) 
{
    std::ifstream inFile(filename);
    if (!inFile) 
    {
        std::cerr << "Error opening file for reading: " << filename << std::endl;
        return {};
    }

    std::vector<LibEntry> entries;
    std::string line;

    while (std::getline(inFile, line)) 
    {
        entries.push_back(LibEntry(line));
    }
    inFile.close();
    return entries;
}


#if false
    /* Deprecated, use force_update_registry followed by filter_registry */
  void update_registry( Delay D, uint32_t func, uint64_t hash,
      std::unordered_map <uint32_t, std::vector<LibEntry>> & functions,
      std::unordered_map<uint64_t, Node> & nodemap)
  {
    bool add = true;
    std::vector<LibEntry> to_delete;

    auto it = functions.find(func);
    if (it != functions.end()) 
    {
      std::unordered_map<Delay<N>, uint64_t>& inner_map = it->second;
      // Iterate over the entries in the inner map
      for (const auto& [other_D, other_hash] : inner_map) 
      {
        if (D.strictly_better_than(other_D))        // new delay pattern dominates the old delay pattern
        {
          to_delete.push_back(other_D); 
        }
        else if (other_D.strictly_better_than(D))   // old delay pattern dominates the new delay pattern
        {
          to_delete.push_back(D);
          add = false;
        }
        else if (D == other_D)                      // the delay patterns are equal, look at the cost
        {
          Node &       node = nodemap[      hash];
          Node & other_node = nodemap[other_hash];
          add = (node.cost < other_node.cost);
        }
      }
    }
    if (add)
    {
      // fmt::print("{}\n", node.cost);
      functions[func][D] = hash;
    }
    for (auto & bad_D : to_delete)
    {
      functions[func].erase(bad_D);
    }
  }
#endif

void force_update_registry( Delay new_D, uint16_t func, uint64_t new_hash, std::vector<uint8_t> pi_order, std::unordered_map<uint16_t, std::unordered_map<Delay, std::vector<LibEntry>>> & functions, std::unordered_map<uint64_t, Node> & nodemap)
{
  bool add = true;
  // Iterate over the entries in the inner map
  for (const auto& [D, entries] : functions[func]) 
  {
    Node & new_node = nodemap[new_hash];
    // New delay pattern dominates the old delay pattern.
    // No processing needed now
    /*
    if (new_D.strictly_better_than(D)){}      
    // Old delay pattern dominates the new delay pattern
    // only add new node if there is no existing node with better cost
    */
    if (D.strictly_better_than(new_D))   
    {
      for (const auto & entry : entries)
      {
        Node & old_node = nodemap[entry.hash];
        if (new_node.cost >= old_node.cost)
        {
          add = false;
          break;
        }
      }
    }      
    // the delay patterns are equal, add for now
    /* else if (new_D == D) {} */
  }
  if (add)
  {
    std::string name = fmt::format("0x{0:04x}_{1}", func, fmt::join(new_D.delays, ""));
    LibEntry new_entry = LibEntry(name, new_D, new_hash, pi_order, nodemap);
    auto shash = new_entry.struct_hash;
    if (std::find(shash.begin(), shash.end(), INF) != shash.end())
    {
      return;
    }
    fmt::print("{0:04x} HASH : new_entry : {1}\n", func, fmt::join(new_entry.struct_hash, ","));
    for (LibEntry entry: functions[func][new_D])
    {
      if (new_entry.struct_hash == entry.struct_hash)
      {
        // fmt::print("Duplicate structure found. Skipping\n");
        return;
      }
    }
    functions[func][new_D].push_back(new_entry); // fmt::print("{}\n", node.cost);
  }
}

void filter_registry( std::unordered_map<uint16_t, std::unordered_map<Delay, std::vector<LibEntry>>> & functions, std::unordered_map<uint64_t, Node> & nodemap)
{

  for (auto & [func, inner_map] : functions)
  {
    fmt::print("Filtering {0:04x}={0:016b}={0}\n", func);
    // First filter each delay pattern based on cost
    for (auto & [D, entries] : inner_map) 
    {
      fmt::print("\tSorting pattern {0} : {1} entries\n", fmt::join(D.delays, ","), entries.size());
      std::sort(
        entries.begin(), 
        entries.end(), 
        [&](const LibEntry& a, const LibEntry& b) 
        {
          return nodemap[a.hash].cost < nodemap[b.hash].cost;
        }
      );
      
      if (entries.size() < MAX_PER_FUNC)
      {
        break;
      }
      // Leave at least <MAX_PER_FUNC> items based on cost 
      auto worst_acceptable_cost = nodemap[entries[MAX_PER_FUNC].hash].cost;

      fmt::print("\tRemoving the entries with cost > {0}\n", worst_acceptable_cost);
      for (auto e_it = entries.begin() + MAX_PER_FUNC; e_it != entries.end(); ++e_it)
      {
          auto current_cost = nodemap[e_it->hash].cost;
          if (current_cost > worst_acceptable_cost)
          {
              // Erase elements from iter onwards
              entries.erase(e_it, entries.end());
              break;
          }
      }
    }
    /* At this point, each delay pattern should have only a few entries */

    // std::unordered_map<Delay, std::vector<LibEntry>>::iterator to_delete

    fmt::print("\tPairwise comparison of delay patterns\n");
    std::vector<Delay> patterns_to_delete;
    for (auto it1 = inner_map.begin(); it1 != inner_map.end(); ++it1) 
    {
      const Delay & Di = it1->first;
      std::vector<LibEntry>& entries_i = it1->second;
      for (auto it2 = std::next(it1); it2 != inner_map.end(); ++it2) 
      {
        const Delay & Dj = it2->first;
        std::vector<LibEntry>& entries_j = it2->second;

        fmt::print("\t{0} vs. {1}\n", fmt::join(Di.delays, ","), fmt::join(Dj.delays, ","));
        if (entries_i.empty()) // entries_i can become empty during the previous iterations
        {
          patterns_to_delete.push_back(Di);
          break;
        }
        if (entries_j.empty()) 
        {
          patterns_to_delete.push_back(Dj);
          continue;
        }
        auto min_cost_i = entries_i.front().cost(nodemap);
        auto min_cost_j = entries_j.front().cost(nodemap);

        if (Di.strictly_better_than(Dj)) 
        /* Only retain those Dj entries whose cost is smaller than smallest-cost entry in Di */
        {
          fmt::print("\t\t{0} is strictly better than {1}\n", fmt::join(Di.delays, ","), fmt::join(Dj.delays, ","));
          for (auto it = entries_j.begin(); it < entries_j.end(); ++it)
          {
            if (min_cost_i <= it->cost(nodemap))
            {
              entries_j.erase(it, entries_j.end());
              break;
            }
          }
        }
        else if (Dj.strictly_better_than(Di)) 
        /* Only retain those Di entries whose cost is smaller than smallest-cost entry in Dj */
        {
          fmt::print("\t\t{1} is strictly better than {0}\n", fmt::join(Di.delays, ","), fmt::join(Dj.delays, ","));
          for (auto it = entries_i.begin(); it < entries_i.end(); ++it)
          {
            if (min_cost_j <= it->cost(nodemap))
            {
              entries_i.erase(it, entries_i.end());
              break;
            }
          }
        }
      }
    }
    for (const Delay & key : patterns_to_delete)
    {
      assert(inner_map[key].empty());
      inner_map.erase(key);
    }
  }
}

// Function to write unordered_map to CSV file
void writeCSV(const std::string& filename, const std::unordered_map<std::string, int>& map) 
{
  std::ofstream outfile(filename);    // Open the output file stream
  outfile << "Key,Value\n";           // Write the header row
  for (const auto& pair : map) 
  {
    outfile << pair.first << "," << pair.second << "\n";    // Write each key-value pair to a new row
  }
  outfile.close();                    // Close the output file stream
}

