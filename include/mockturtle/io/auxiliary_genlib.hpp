#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <fmt/format.h>
#include <mockturtle/algorithms/nodes.hpp>

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
  std::string chars;  
  /*
    PI order is [1,0,3] means that inputs [a,b,c] to a node in kLUT ntk correspond to PIs [0x3333, 0x5555, 0x00FF], respectively. 
  */
  std::vector<UI> struct_hash;

  // Copy constructor
  LibEntry(const LibEntry& other) :
    name(other.name),
    delays(other.delays),
    hash(other.hash),
    chars(other.chars), 
    struct_hash(other.struct_hash) {}

  // Explicit constructor
  explicit LibEntry
  (
    const std::string _name,
    const Delay& _delays, 
    const uint64_t _hash, 
    const std::string& _chars,
    phmap::flat_hash_map<uint64_t, Node> & nodemap
  ) :
    name(_name),
    delays(_delays),
    hash(_hash),
    chars(_chars) 
    {
      struct_hash = nodemap[hash].structural_hash(nodemap);
    }

  explicit LibEntry(const std::string& csv_line) 
  {
    std::istringstream iss(csv_line);
    std::string token;

    // fmt::print("Reading the line {0}\n", csv_line);

    std::getline(iss, name, ',');
    // fmt::print("\tname: {0}\n", name);

    std::getline(iss, token, ',');
    std::stringstream ss(token);
    std::string token_delay;

    while (std::getline(ss, token_delay, '|')) {
        delays.delays.push_back(std::stoi(token_delay));
    }

    // while (std::getline(iss, token, '|')) 
    // {
    //     delays.delays.push_back(static_cast<uint8_t>(std::stoi(token)));
    // }
    // // read last delay
    // std::getline(iss, token, ',');
    // delays.delays.push_back(static_cast<uint8_t>(std::stoi(token)));

    // fmt::print("\tdelays: ({0})\n", fmt::join(delays.delays, ","));

    std::getline(iss, token, ',');
    // fmt::print("\thash: {0}\n", token);
    hash = std::stoull(token);
    std::getline(iss, chars, ',');
    // fmt::print("\tchars: {0}\n", chars);
  }

  // Assignment operator
  LibEntry& operator=(const LibEntry& other) 
  {
    name = other.name;
    delays = other.delays;
    hash = other.hash;
    chars = other.chars;
    struct_hash = other.struct_hash;
    return *this;
  };
  UI cost(phmap::flat_hash_map<uint64_t, Node> & nodemap) const
  {
    return nodemap[hash].cost;
  }
  // UI cost(phmap::flat_hash_map<uint64_t, Node> & nodemap) const
  // {
  //   return nodemap[hash].cost;
  // }

  std::string write() const
  {
    return fmt::format(
      "{0}\n",
      to_str()
    );
  }
  std::string to_str() const
  {
    return fmt::format(
      "{0},{1},{2},{3}",
      name,
      fmt::join(delays.delays, "|"),
      hash,
      chars
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

phmap::flat_hash_map<std::string, LibEntry> read_LibEntry_map(const std::string& filename)
{
  std::ifstream inFile(filename);
  if (!inFile) 
  {
    std::cerr << "Error opening file for reading: " << filename << std::endl;
    return {};
  }

  phmap::flat_hash_map<std::string, LibEntry> entries;
  std::string line;
  while (std::getline(inFile, line)) 
  {
    // entries.push_back(LibEntry(line));
    LibEntry entry = LibEntry(line);
    entries.emplace(entry.name, entry);
  }
  inFile.close();
  return entries;
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
