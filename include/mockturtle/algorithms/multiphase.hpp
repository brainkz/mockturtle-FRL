
#pragma once

#include <mockturtle/io/auxiliary_genlib.hpp>
#include <mockturtle/algorithms/nodes.hpp>
#include <mockturtle/utils/misc.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/rsfq_view.hpp>
#include <mockturtle/views/fanout_view.hpp>


typedef uint32_t glob_phase_t;

typedef mockturtle::klut_network klut;
typedef mockturtle::xag_network   xag;
typedef mockturtle::xmg_network   xmg;
typedef mockturtle::mig_network   mig;
typedef mockturtle::aig_network   aig;
typedef uint64_t node_t;

enum GateType : uint8_t 
{
    PI_GATE = 0u,
    AA_GATE = 1u,
    AS_GATE = 2u,
    SA_GATE = 3u,
    T1_GATE = 4u
};

const std::vector<std::string> GATE_TYPE { "PI", "AA", "AS", "SA", "T1" }; //, "PO"

union NodeData 
{
  struct 
  {
    unsigned int sigma : 29;
    unsigned int type : 3;
  };
  uint32_t value;
  NodeData() : value(0) {}
  NodeData(uint32_t _sigma, uint8_t _type) : sigma(_sigma), type(_type) {}
  NodeData(uint32_t _value) : value(_value) {}
};

struct Path
{
  std::set<klut::signal> sources;   // AS/SA gates
  std::set<klut::signal> internals; // AA    gates
  std::set<klut::signal> targets;   // AS/SA gates
  Path(const std::set<klut::signal> & _sources, const std::set<klut::signal>& _internals, const std::set<klut::signal>& _targets)
    : sources(_sources), internals(_internals), targets(_targets) {}
  Path() : sources({}), internals({}), targets({}) {}
  
  void absorb(Path & other)
  {
    sources.insert(other.sources.begin(), other.sources.end());
    internals.insert(other.internals.begin(), other.internals.end());
    targets.insert(other.targets.begin(), other.targets.end());
  }

  void print() const
  {
    fmt::print("Path from [{}]\n\tvia [{}]\n\tto [{}]\n", fmt::join(sources, ","), fmt::join(internals, ","), fmt::join(targets, ","));
  }
  std::string format() const
  {
    return fmt::format("Path from [{}]\n\tvia [{}]\n\tto [{}]\n", fmt::join(sources, ","), fmt::join(internals, ","), fmt::join(targets, ","));
  }

  void print_bfs(const klut & ntk) const
  {
    fmt::print("Path from [{}]\n\tvia [{}]\n\tto [{}]\n", fmt::join(sources, ","), fmt::join(internals, ","), fmt::join(targets, ","));
    std::deque<klut::signal> queue { targets.begin(), targets.end() };
    while (!queue.empty())
    {
      const klut::signal & node = queue.front();
      queue.pop_front();

      const std::string kind { _kind(node) } ;

      const std::vector<klut::signal> fanins { preds(node, ntk) };

      const NodeData node_data { ntk.value(node) };

      fmt::print("\t[{} {} {}]\tσ={}\tfanins : {{ {} }} \n", GATE_TYPE[node_data.type], kind, node, static_cast<int>(node_data.sigma), fmt::join(fanins, ", "));

      for (const klut::signal & fanin : fanins)
      {
        queue.push_back(fanin);
      }
    }
    return;
  }

  std::vector<klut::signal> preds(const klut::signal & sig, const klut & ntk) const
  {
    if ( sources.count(sig) != 0)
    {
      return {};
    }

    std::vector<klut::signal> predecessors;

    ntk.foreach_valid_fanin(sig, [&](const klut::signal & parent)
    {
      if ( internals.count(parent) != 0 || sources.count(parent) != 0 )
      {
        predecessors.push_back( parent );
      }
    });
    return predecessors;
  }

  std::vector<klut::signal> src_int() const
  {
    std::vector<klut::signal> out;
    out.insert(out.end(), sources.begin(), sources.end());
    out.insert(out.end(), internals.begin(), internals.end());
    return out;
  }
  std::vector<klut::signal> int_tgt() const
  {
    std::vector<klut::signal> out;
    out.insert(out.end(), internals.begin(), internals.end());
    out.insert(out.end(), targets.begin(), targets.end());
    return out;
  }

  std::string _kind( const klut::signal & node ) const
  {
    if ( targets.count(node) > 0 )  
    { 
      return "Target";  
    }
    else if ( internals.count(node) > 0 )  
    { 
      return "Internal";  
    }
    else if ( sources.count(node) > 0 )  
    { 
      return "Source";  
    }
    else 
    {
      throw;
    }
  }

  /// @brief Find paths of signals in Path object from target signals to their sources.
  /// This function performs a depth-first search starting from target signals and follows
  /// their predecessors until no more predecessors are found, creating threads of signals.
  /// 
  /// @param ntk The logic network.
  /// @return A vector of vectors, where each inner vector represents a thread of signals
  ///         from a target signal to its sources.
  std::vector<std::vector<klut::signal>> path_threads(const klut & ntk, bool verbose = false) const
  {
    // Initialize two vectors to store incomplete and completed threads
    std::vector<std::vector<klut::signal>> incomplete_threads;
    std::vector<std::vector<klut::signal>> done_threads;

    // Start a thread for each target signal
    for (const klut::signal & tgt : targets)
    {
      incomplete_threads.push_back({ tgt });
    }

    DEBUG_PRINT("[path_threads] CREATED {} incomplete_threads\n", incomplete_threads.size());
    for (auto thread: incomplete_threads)
    {
      DEBUG_PRINT("[path_threads]\tincomplete thread: [{}]\n", fmt::join(thread, ","));
    }

    // Continue until there are no incomplete threads left
    while (!incomplete_threads.empty())
    {
      // Take the last thread from the incomplete_threads vector
      std::vector<klut::signal> thread = incomplete_threads.back();
      incomplete_threads.pop_back();

      DEBUG_PRINT("[path_threads] processing thread: [{}]\n", fmt::join(thread, ","));

      // Get the gate signal at the end of the current thread
      const klut::signal & gate = thread.back();

      // Find predecessors of the gate signal in the network
      std::vector<klut::signal> predecessors { preds(gate, ntk) };

      DEBUG_PRINT("[path_threads]\tpredecessors: [{}]\n", fmt::join(predecessors, ","));

      // If there are no predecessors, this thread is complete
      if (predecessors.empty())
      {
        done_threads.push_back( thread );
        for (auto thread: done_threads)
        {
          DEBUG_PRINT("[path_threads] thread complete: [{}]\n", fmt::join(thread, ","));
        }
      }
      else
      {
        // If there are predecessors, extend the thread with each predecessor
        for ( const klut::signal & predecessor : predecessors )
        {
          incomplete_threads.push_back( thread );
          incomplete_threads.back().push_back(predecessor);
          DEBUG_PRINT("[path_threads]\tnew thread: [{}]\n", fmt::join(incomplete_threads.back(), ","));
        }
      }
    }

    // Return the completed threads
    return done_threads;
  } 
};


std::tuple<klut, int64_t> decompose_to_klut(mockturtle::binding_view<klut> src, phmap::flat_hash_map<ULL, Node> nodemap, phmap::flat_hash_map<std::string, LibEntry> entries, const std::array<int, 12> COSTS_MAP, bool verbose = false)
{
  phmap::flat_hash_map<klut::signal, klut::signal> src2tgt;
  int64_t area = 0;
  
  klut tgt;
  src.foreach_pi( [&]( const klut::signal & src_pi ) 
  {
    klut::signal tgt_pi = tgt.create_pi();
    tgt.set_value(tgt_pi, NodeData(0, PI_GATE).value);
    src2tgt.emplace(src_pi, tgt_pi);
  } );

  // The bindings are replaced in forward topological order.
  std::vector<klut::signal> src_po;
  src_po.reserve(src.num_pos());
  src.foreach_po( [&] (auto src_n)
  {
    src_po.push_back(src_n);
  } );

  std::unordered_set<klut::signal> OR_replacement_candidates;
  std::vector<klut::signal> XOR_gates;

  auto num_node_src = src.size();
  auto ctr = 0u;
  mockturtle::topo_view<klut>( src ).foreach_node( [&]( auto src_node ) 
  {
    DEBUG_PRINT("Processing node {0} ({1} out of {2})\r", src_node, ++ctr, num_node_src);
    if (ctr == num_node_src)
    {
      DEBUG_PRINT("\n");
    }
    if ( !src.has_binding( src_node ) )
    {
      return;
    }

    auto const& g = src.get_binding( src_node );

    // handing constants
    if (g.name == "one" && g.expression == "CONST1")
    {
      klut::signal tgt_sig = tgt.get_constant( true );
      src2tgt.emplace(src_node, tgt_sig);
      return;
    }
    else if (g.name == "zero" && g.expression == "CONST0")
    {
      klut::signal tgt_sig = tgt.get_constant( false );
      src2tgt.emplace(src_node, tgt_sig);
      return;
    }
    
    LibEntry entry = entries.at(g.name);
    Node & root = nodemap[entry.hash];

    // Get the topological order of internal nodes (i.e., nodes inside the cell)
    std::vector<ULL> topo_order;
    root.topo_sort(nodemap, topo_order, true);

    std::vector<klut::signal> tgt_fanins;
    src.foreach_fanin(src_node, [&](const auto & src_fanin)
    {
      tgt_fanins.push_back(src2tgt.at(src_fanin));
      // fmt::print("\tRecorded fanin of src_node {0}:\n\t\tsrc_fanin: {1}\n\t\ttgt_fanin: {2}\n", src_node, src_fanin, tgt_fanins.back().data);
    } );

    phmap::flat_hash_map<ULL, klut::signal> node2tgt;

    for (ULL hash : topo_order)
    {
      Node & node = nodemap.at(hash);
      if (node.last_func == fPI)
      {
        // determine which index of the PI is needed 
        char letter = node.pi_letter();

        auto pi_idx = std::find(entry.chars.begin(), entry.chars.end(), letter) - entry.chars.begin();
        node2tgt.emplace(hash, tgt_fanins[pi_idx]);
      }
      else if (node.last_func == fDFF)
      {
        klut::signal parent_tgt = node2tgt.at(node.parent_hashes.front());
        node2tgt.emplace(hash, parent_tgt);
        //  nothing else to do, no new constraints, no contribution to objective function
      }
      else if (node.last_func == fNOT)
      {
        klut::signal parent_tgt = node2tgt.at(node.parent_hashes.front());
        klut::signal tgt_sig = tgt.create_not( parent_tgt );
        tgt.set_value(tgt_sig, NodeData(0, AS_GATE).value);
        node2tgt.emplace(hash, tgt_sig);
        area += COSTS_MAP[fNOT];
        DEBUG_PRINT("ADDED NOT = {}\n", area);
      }
      else if (node.last_func == fAND)
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_and(parent_tgt_1, parent_tgt_2);
        tgt.set_value(tgt_sig, NodeData(0, SA_GATE).value);
        // fmt::print("Created node n{0} = AND({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        area += COSTS_MAP[fAND];
        DEBUG_PRINT("ADDED AND = {}\n", area);
      }
      else if (node.last_func == fOR) 
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_or(parent_tgt_1, parent_tgt_2);
        tgt.set_value(tgt_sig, NodeData(0, SA_GATE).value);
        // fmt::print("Created node n{0} = OR({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        OR_replacement_candidates.emplace(tgt_sig);
        area += COSTS_MAP[fOR];
        DEBUG_PRINT("ADDED OR  = {}\n", area);
      }
      else if (node.last_func == fCB)
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_or(parent_tgt_1, parent_tgt_2);
        tgt.set_value(tgt_sig, NodeData(0, AA_GATE).value);
        // fmt::print("Created node n{0} = CB({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        area += COSTS_MAP[fCB];
        DEBUG_PRINT("ADDED CB  = {}\n", area);
      }
      else if (node.last_func == fXOR)
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_xor(parent_tgt_1, parent_tgt_2);
        tgt.set_value(tgt_sig, NodeData(0, AS_GATE).value);
        // fmt::print("Created node n{0} = XOR({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        XOR_gates.push_back(tgt_sig);
        area += COSTS_MAP[fXOR];
        DEBUG_PRINT("ADDED XOR = {}\n", area);
      }
      // else if 
      else
      {
        throw "Unsupported function";
      }
    }
    ULL root_hash = topo_order.back();
    klut::signal tgt_sig = node2tgt.at(root_hash);
    src2tgt.emplace(src_node, tgt_sig);
  } );

  for (klut::signal const & xor_sig : XOR_gates)
  {
    tgt.foreach_fanin(xor_sig, [&] (const klut::signal & fanin)
    {
      // if the OR gate is before XOR, do not replace it with the CB
      auto it = std::find(OR_replacement_candidates.begin(), OR_replacement_candidates.end(), fanin);
      if (it != OR_replacement_candidates.end())
      {
        OR_replacement_candidates.erase(it);
      }
    });
  }

  for (klut::signal const & sig : OR_replacement_candidates)
  {
    tgt.set_value(sig, NodeData(0, AA_GATE).value);
    area += COSTS_MAP[fCB];
    area -= COSTS_MAP[fOR];
    DEBUG_PRINT("REPLACED OR WITH CB = {}\n", area);
  }

  src.foreach_po([&](auto const & src_po)
  {
    tgt.create_po(src2tgt.at(src_po));
  } );
  
  return std::make_tuple(tgt, area);
}


template <typename Ntk>
glob_phase_t latest_fanin_phase(const Ntk & ntk, const typename Ntk::signal & node, const uint8_t n_phases, const uint8_t type, const bool verbose = false)
{
  bool valid = false;
  uint32_t phase = 0u;
  ntk.foreach_fanin(node, [&] ( const typename Ntk::signal & parent )
  {
    if ( ntk.is_constant( parent ) )
    {
      DEBUG_PRINT("\t\tfanin {} is constant, skipping...\n", parent);
      return;
    }
    valid = true;
    NodeData parent_data = ntk.value( parent );

    // if d==false, SA gate can be directly connected to fanin
    //cannot directly connect PI (convention)
    //cannot directly connect split signal 
    //cannot directly connect SA/AA gates to SA gates
    int d = (type == SA_GATE) && ( (ntk.is_pi(parent)) || (ntk.fanout_size(parent) > 1) || (parent_data.type != AS_GATE) );

    phase = generic_max( phase, parent_data.sigma + d );

    if (verbose)
    {
      unsigned int gp = static_cast<unsigned int>(parent_data.sigma);
      DEBUG_PRINT("\t\tfanin {} ɸ={} [S={}, φ={}]\n", parent, gp, gp/n_phases, gp%n_phases);
    }
  });
  assert(valid);

  // DEBUG_PRINT("\t{} GATE {} placed at ɸ={} [S={}, φ={}]\n",  node, phase, phase/n_phases, phase%n_phases);

  return phase;
}

/// @brief Assigns stages to nodes based on stage assignment. If the assignment is not found, assigns a stage greedily
/// @param ntk 
/// @param n_phases 
/// @param phase_assignment 
/// @param verbose 
void assign_sigma(const klut & ntk, const phmap::flat_hash_map<unsigned int, unsigned int> & phase_assignment, const bool verbose = false)
{
  mockturtle::topo_view<klut> ntk_topo ( ntk );

  ntk_topo.foreach_node([&] ( const klut::signal & node ) 
  {
    if ( ntk_topo.is_constant( node ) )
    {
      return;
    }

    if (phase_assignment.count(node) != 0) // there is a precalculated phase assignment
    {
      uint8_t node_type = NodeData( ntk.value( node ) ).type;

      NodeData node_data;
      node_data.type = node_type;
      node_data.sigma = phase_assignment.at( node );

      ntk.set_value(node, node_data.value);
      DEBUG_PRINT("{} GATE {}, σ={}\n", GATE_TYPE.at(node_type), node, static_cast<int>(node_data.sigma));

      ntk_topo.foreach_valid_fanin(node, [&](const klut::signal & fi_node)
      {
        NodeData fi_node_data;
        DEBUG_PRINT("\tFANIN {}, σ={}\n", fi_node, static_cast<int>(fi_node_data.sigma));
        // node_data.sigma = phase_assignment.at( node );
      });
    }
  });
}


/// @brief Assigns stages to nodes based on stage assignment. If the assignment is not found, assigns a stage greedily
/// @param ntk 
/// @param n_phases 
/// @param phase_assignment 
/// @param verbose 
void greedy_ntk_assign(const klut & ntk, const uint8_t n_phases, const phmap::flat_hash_map<unsigned int, unsigned int> & phase_assignment, const bool verbose = false)
{
  mockturtle::topo_view<klut> ntk_topo ( ntk );

  ntk_topo.foreach_node([&] ( const klut::signal & node ) 
  {
    if ( ntk_topo.is_constant( node ) )
    {
      return;
    }

    if ( ntk_topo.is_pi( node ) )
    {
      auto ct = phase_assignment.count(node);
      uint32_t sigma = (ct != 0) ? phase_assignment.at( node ) : 0;
      NodeData node_data {sigma, AS_GATE};
      ntk.set_value(node, node_data.value);

      DEBUG_PRINT("PI {} placed at ɸ=0 [S=0, φ=0]\n", node);
      return;
    }

    // if (verbose) fmt::print("{} GATE {}:\n", GATE_TYPE.at(node_params.type), node);

    uint8_t node_type = NodeData( ntk.value( node ) ).type;
    DEBUG_PRINT("{} GATE {}:\n", GATE_TYPE.at(node_type), node);

    NodeData node_data;
    node_data.type = node_type;

    auto ct = phase_assignment.count(node);
    if (ct != 0) // there is a precalculated phase assignment
    {
      node_data.sigma = phase_assignment.at( node );
      ntk.set_value(node, node_data.value);
    }
    else
    {
      if ( node_type == AA_GATE )
      {
        node_data.sigma = latest_fanin_phase(ntk_topo, node, n_phases, AA_GATE, verbose);
        uint32_t sigma = static_cast<uint32_t>(node_data.sigma);
        DEBUG_PRINT("\tAA GATE {} placed at ɸ={} [S={}, φ={}]\n", node, sigma, sigma/n_phases, sigma%n_phases);
      }
      else if ( node_type == SA_GATE )
      {
        node_data.sigma = latest_fanin_phase(ntk_topo, node, n_phases, SA_GATE, verbose);
        uint32_t sigma = static_cast<uint32_t>(node_data.sigma);
        DEBUG_PRINT("\tSA GATE {} placed at ɸ={} [S={}, φ={}]\n", node, sigma, sigma/n_phases, sigma%n_phases);
      }
      else if ( node_type == AS_GATE )
      {
        node_data.sigma = latest_fanin_phase(ntk_topo, node, n_phases, SA_GATE, verbose) + 1;
        uint32_t sigma = static_cast<uint32_t>(node_data.sigma);
        DEBUG_PRINT("\tAS GATE {} placed at ɸ={} [S={}, φ={}]\n", node, sigma, sigma/n_phases, sigma%n_phases);
      }
      else 
      {
        std::vector<klut::signal> fanins;
        ntk.foreach_fanin(node, [&](const klut::signal & fanin){ fanins.push_back(fanin); });
        fmt::print("\t GATE {0} : fanins[{1}], type[{2}], func[{3}=0x{3:x}=0b{3:b}]\n", node, fmt::join(fanins, ","), node_type, ntk.node_function(node)._bits[0]);
        throw "Unsupported gate type!";
      }
    }
  });
}

/// @brief compares stages of two klut nodes. 
/// @param a - first klut signal to be compared
/// @param b - first klut signal to be compared
/// @param ntk - klut network
/// @return 
bool phase_ntk_comparison(const klut::signal & a, const klut::signal & b, const klut & ntk )
{
  if (ntk.is_constant(a))
  {
    return true;
  }
  else if (ntk.is_constant(b))
  {
    return false;
  }
  NodeData a_data = ntk.value(a);
  NodeData b_data = ntk.value(b);

  return a_data.sigma < b_data.sigma;
}


// Function to insert splitter nodes in a KLUT network.
void splitter_ntk_insertion(klut & ntk, const bool verbose = false)
{
  // Lambda function for comparing the phases of two signals in the network.
  auto phase_comp = [&](const klut::signal & a, const klut::signal & b)
  {
    return phase_ntk_comparison(a, b, ntk);
  };

  // Create a view of the network that provides access to fanout information.
  auto ntk_fo = mockturtle::fanout_view<klut>(ntk);

  // For each node in the fanout view:
  ntk_fo.foreach_node([&](const klut::signal & node)
  {
    if ( ntk_fo.is_dangling( node ) )
    {
      return;
    }

    // Get the number of fanouts for the current node.
    uint32_t fo_size{ 0u };
    ntk_fo.foreach_fanout( node, [&]( auto const& no ) {
      if ( !ntk_fo.is_dangling( no ) )
      {
        ++fo_size;
      }
    } );
    ntk._storage->nodes[node].data[0].h1 = fo_size;

    DEBUG_PRINT("\t[NODE {}] FANOUT SIZE = {}\n", node, fo_size);
    // If the current node is a constant or it has fanout ≤ 1, skip to the next node.
    if (ntk_fo.is_constant(node) || fo_size <= 1)
    {
      return;
    }

    // Populate the fanouts vector.
    std::vector<klut::signal> fanouts;
    fanouts.reserve(fo_size);
    ntk_fo.foreach_fanout(node, [&](const klut::signal & fo_node)
    {
      if ( ntk_fo.is_dangling( fo_node ) )
      {
        return;
      }

      fanouts.push_back(fo_node);
      DEBUG_PRINT("\t\t[NODE {}] ADDING FANOUT\n", node, fo_node);
    });

    // Fix the fanout count (bugged fanouts_size()?)
    if ( fanouts.size() != fo_size )
    {
      ntk._storage->nodes[node].data[0].h1 = fanouts.size();
    }

    // Sort the fanouts using the phase comparison function.
    std::sort(fanouts.begin(), fanouts.end(), phase_comp);
    DEBUG_PRINT("\t[NODE {}] SORTED FANOUTS:\n", node);
    printVector(fanouts, 2);

    // Create [fo_size - 1] splitter nodes.
    klut::signal last_spl = node;
    std::vector<klut::signal> splitters;
    splitters.reserve( fanouts.size() - 1 );
    for (auto it = fanouts.begin(); it < fanouts.end() - 1; it++)
    {
      DEBUG_PRINT("\t\t[NODE {}] LAST SPL: {}\n", node, last_spl);
      // Copy sigma and type data from the first fanout for the splitter node data.
      NodeData spl_data = ntk.value(*it);
      // Change gate type to AA for the splitter node.
      spl_data.type = AA_GATE;
      // Create a new splitter node connected to 'last_spl'.

      DEBUG_PRINT("\t\t[NODE {}] CREATING SPL FOR {}\n", node, *it);
      DEBUG_PRINT("\t\t[NODE {}] LAST_SPL {} FANOUT BEFORE: {}\n", node, last_spl, ntk.fanout_size(last_spl));
      const klut::signal spl = ntk._create_node({last_spl}, 2, last_spl);
      ntk.set_value(spl, spl_data.value);
      DEBUG_PRINT("\t\t[NODE {}] CREATED SPL {}\n", node, spl);
      DEBUG_PRINT("\t\t[NODE {}] LAST_SPL {} FANOUT AFTER: {}\n", node, last_spl, ntk.fanout_size(last_spl));
      DEBUG_PRINT("\t\t[NODE {}] SPL {} FANIN: {}\n", node, spl, ntk.fanin_size(spl));

      // auto n = ntk._storage->nodes[*it];
      // auto& preds = n.children;

      // Update the connections to reflect the splitter's presence.
      DEBUG_PRINT("\t\t[NODE {}, LAST_SPL {}] UPDATING CONNECTIONS\n", node, last_spl);
      // for (auto& pred : preds)
      for (auto pred_it = ntk._storage->nodes[*it].children.begin();
                pred_it < ntk._storage->nodes[*it].children.end(); pred_it++ )
      {
        if ( ntk.is_dangling( pred_it->data ) )
        {
          continue;
        }

        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PRED {}\n", node, last_spl, pred_it->data);
        if (pred_it->data == node)
        {
          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] RECORDING OLD PREDS\n", node, last_spl);
          // Store the previous connections.
          std::vector<klut::signal> old_preds(ntk._storage->nodes[*it].children.size());
          std::transform(ntk._storage->nodes[*it].children.begin(), ntk._storage->nodes[*it].children.end(), old_preds.begin(), [](auto c) { return c.index; });
          printVector(old_preds, 4);

          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS BEFORE\n", node, last_spl);
          for (const auto& entry : ntk._storage->nodes[*it].children)  { DEBUG_PRINT("\t\t\t\t{}\n", entry.data); }
          pred_it->data = spl;                             // Replace the connection with the splitter.
          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS AFTER\n", node, last_spl);
          for (const auto& entry : ntk._storage->nodes[*it].children)  { DEBUG_PRINT("\t\t\t\t{}\n", entry.data); }

          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] SPL {} FANOUT BEFORE: {}\n", node, last_spl, spl, static_cast<int>(ntk._storage->nodes[spl].data[0].h1));
          DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(spl));
          ntk._storage->nodes[spl].data[0].h1++;  // Increment fan-out of the splitter.
          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] SPL {} FANOUT  AFTER: {}\n", node, last_spl, spl, static_cast<int>(ntk._storage->nodes[spl].data[0].h1));
          DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(spl));

          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] NODE {} FANOUT BEFORE: {}\n", node, last_spl, node, static_cast<int>(ntk._storage->nodes[node].data[0].h1));
          DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(node));
          ntk._storage->nodes[node].data[0].h1--; // Decrement fan-out of the current node.
          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] NODE {} FANOUT  AFTER: {}\n", node, last_spl, node, static_cast<int>(ntk._storage->nodes[node].data[0].h1));
          DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(node));

          // // Notify listeners of the modification.
          // for (auto const& fn : ntk._events->on_modified)
          // {
          //   (*fn)(*it, old_preds);
          // }
          break;
        }
      }
      last_spl = spl;
      splitters.push_back(spl);
    }

    // Process the last fanout.
    // auto& preds = ntk._storage->nodes[fanouts.back()].children;
    DEBUG_PRINT("\t\t[NODE {}, LAST_SPL {}] UPDATING CONNECTIONS TO LAST FANOUT {}\n", node, last_spl, fanouts.back());
    // for (auto& pred : preds)
    for (auto pred_it = ntk._storage->nodes[fanouts.back()].children.begin();
              pred_it < ntk._storage->nodes[fanouts.back()].children.end(); pred_it++ )
    {
      if ( ntk.is_dangling( pred_it->data ) )
      {
        continue;
      }

      if (pred_it->data == node)
      {
        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] RECORDING OLD PREDS\n", node, last_spl);
        // Store the previous connections.
        std::vector<klut::signal> old_preds(ntk._storage->nodes[fanouts.back()].children.size());
        std::transform(ntk._storage->nodes[fanouts.back()].children.begin(), ntk._storage->nodes[fanouts.back()].children.end(), old_preds.begin(), [](auto c) { return c.index; });
        printVector(old_preds, 4);

        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS BEFORE\n", node, last_spl);
        for (const auto& entry : ntk._storage->nodes[fanouts.back()].children)  { DEBUG_PRINT("\t\t\t\t{}\n", entry.data); }
        pred_it->data = last_spl;                            // Replace the connection with the last splitter.
        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS AFTER\n", node, last_spl);
        for (const auto& entry : ntk._storage->nodes[fanouts.back()].children)  { DEBUG_PRINT("\t\t\t\t{}\n", entry.data); }     


        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] SPL {} FANOUT BEFORE: {}\n", node, last_spl, last_spl, static_cast<int>(ntk._storage->nodes[last_spl].data[0].h1));
        DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(last_spl));
        ntk._storage->nodes[last_spl].data[0].h1++; // Increment fan-out of the last splitter.
        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] SPL {} FANOUT  AFTER: {}\n", node, last_spl, last_spl, static_cast<int>(ntk._storage->nodes[last_spl].data[0].h1));
        DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(last_spl));

        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] NODE {} FANOUT BEFORE: {}\n", node, last_spl, node, static_cast<int>(ntk._storage->nodes[node].data[0].h1));
        DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(node));
        ntk._storage->nodes[node].data[0].h1--;     // Decrement fan-out of the current node.
        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] NODE {} FANOUT  AFTER: {}\n", node, last_spl, node, static_cast<int>(ntk._storage->nodes[node].data[0].h1));
        DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(node));

        // Notify listeners of the modification.
        for (auto const& fn : ntk._events->on_modified)
        {
          (*fn)(fanouts.back(), old_preds);
        }
      }
    }

    auto fo_ntk { mockturtle::fanout_view( ntk ) };
    for (const auto & node : fanouts)
    {
      DEBUG_PRINT("\t\t\t node: {}\n", node);
      
      fo_ntk.foreach_fanin( node, [&](const klut::signal & fi_node)
      {
        if ( fo_ntk.is_dangling( fi_node ) )
        {
          return;
        }

        DEBUG_PRINT("\t\t\t\t fanin : {}\n", fi_node);
      });
      fo_ntk.foreach_fanout( node, [&](const klut::signal & fo_node)
      {
        if ( fo_ntk.is_dangling( fo_node ) )
        {
          return;
        }

        DEBUG_PRINT("\t\t\t\t fanout: {}\n", fo_node);
      });
    }
    for (const auto & node : splitters)
    {
      DEBUG_PRINT("\t\t\t spl : {}\n", node);
      
      fo_ntk.foreach_fanin( node, [&](const klut::signal & fi_node)
      {
        if ( fo_ntk.is_dangling( fi_node ) )
        {
          return;
        }

        DEBUG_PRINT("\t\t\t\t fanin : {}\n", fi_node);
      });
      fo_ntk.foreach_fanout( node, [&](const klut::signal & fo_node)
      {
        if ( fo_ntk.is_dangling( fo_node ) )
        {
          return;
        }

        DEBUG_PRINT("\t\t\t\t fanout: {}\n", fo_node);
      });
    }
    // Ensure that the current node's fan-out count is now 1 (since all other fanouts have been replaced by splitters).
    assert(ntk._storage->nodes[node].data[0].h1 == 1);
  });
}

void splitter_ntk_insertion_t1( klut& ntk, phmap::flat_hash_map<klut::signal, klut::signal> const& representatives, const bool verbose = false )
{
  // Lambda function for comparing the phases of two signals in the network.
  auto phase_comp = [&](const klut::signal & a, const klut::signal & b)
  {
    return phase_ntk_comparison(a, b, ntk);
  };

  // Create a view of the network that provides access to fanout information.
  auto ntk_fo = mockturtle::fanout_view<klut>(ntk);

  // For each node in the fanout view:
  ntk_fo.foreach_node( [&](const klut::signal & node) {
    if ( ntk_fo.is_dangling( node ) )
    {
      return;
    }
    
    // Populate the fanouts vector.
    std::vector<klut::signal> fanouts;
    ntk_fo.foreach_fanout( node, [&]( const klut::signal & fo_node ) {
      if ( ntk_fo.is_dangling( fo_node ) )
      {
        return;
      }

      if ( static_cast<NodeData>( ntk.value( fo_node ) ).type == T1_GATE ) 
      {
        /* the current node is an output of a T1 gate */
        if ( representatives.at( fo_node ) != fo_node )
        {
          /* the current node is not a representitive */
          return;
        }
      }

      fanouts.push_back( fo_node );
      DEBUG_PRINT("\t\t[NODE {}] ADDING FANOUT {}\n", node, fo_node);
    } );

    uint32_t fo_size{ fanouts.size() };
    ntk._storage->nodes[node].data[0].h1 = fo_size;
    DEBUG_PRINT("\t[NODE {}] FANOUT SIZE = {}\n", node, fo_size );

    // If the current node is a constant or it has fanout ≤ 1, skip to the next node.
    if (ntk_fo.is_constant(node) || fo_size <= 1)
    {
      return;
    }
    // Sort the fanouts using the phase comparison function.
    std::sort( fanouts.begin(), fanouts.end(), phase_comp );
    DEBUG_PRINT( "\t[NODE {}] SORTED FANOUTS:\n", node );
    printVector( fanouts, 2 );

    // Create [fo_size - 1] splitter nodes.
    klut::signal last_spl = node;
    std::vector<klut::signal> splitters;
    splitters.reserve( fanouts.size() - 1 );
    for ( auto it = fanouts.begin(); it < fanouts.end() - 1; it++ )
    {
      DEBUG_PRINT("\t\t[NODE {}] LAST SPL: {}\n", node, last_spl);
      // Copy sigma and type data from the first fanout for the splitter node data.
      NodeData spl_data = ntk.value( *it );
      // Change gate type to AA for the splitter node.
      spl_data.type = AA_GATE;
      // Create a new splitter node connected to 'last_spl'.

      DEBUG_PRINT("\t\t[NODE {}] CREATING SPL FOR {}\n", node, *it);
      DEBUG_PRINT("\t\t[NODE {}] LAST_SPL {} FANOUT BEFORE: {}\n", node, last_spl, ntk.fanout_size( last_spl ));
      const klut::signal spl = ntk._create_node({last_spl}, 2, last_spl);
      ntk.set_value(spl, spl_data.value);
      DEBUG_PRINT("\t\t[NODE {}] CREATED SPL {}\n", node, spl);
      DEBUG_PRINT("\t\t[NODE {}] LAST_SPL {} FANOUT AFTER: {}\n", node, last_spl, ntk.fanout_size( last_spl ));
      DEBUG_PRINT("\t\t[NODE {}] SPL {} FANIN: {}\n", node, spl, ntk.fanin_size( spl ));

      // auto n = ntk._storage->nodes[*it];
      // auto& preds = n.children;

      // Update the connections to reflect the splitter's presence.
      DEBUG_PRINT("\t\t[NODE {}, LAST_SPL {}] UPDATING CONNECTIONS\n", node, last_spl);
      // for (auto& pred : preds)
      for (auto pred_it = ntk._storage->nodes[*it].children.begin();
                pred_it < ntk._storage->nodes[*it].children.end(); pred_it++ )
      {
        if ( ntk.is_dangling( pred_it->data ) )
        {
          continue;
        }

        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PRED {}\n", node, last_spl, pred_it->data);
        if ( pred_it->data == node )
        {
          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] RECORDING OLD PREDS\n", node, last_spl);
          // Store the previous connections.
          std::vector<klut::signal> old_preds(ntk._storage->nodes[*it].children.size());
          std::transform(ntk._storage->nodes[*it].children.begin(), ntk._storage->nodes[*it].children.end(), old_preds.begin(), [](auto c) { return c.index; });
          printVector(old_preds, 4);

          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS BEFORE\n", node, last_spl);
          for (const auto& entry : ntk._storage->nodes[*it].children)  { DEBUG_PRINT("\t\t\t\t{}\n", entry.data); }
          pred_it->data = spl;                             // Replace the connection with the splitter.
          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS AFTER\n", node, last_spl);
          for (const auto& entry : ntk._storage->nodes[*it].children)  { DEBUG_PRINT("\t\t\t\t{}\n", entry.data); }

          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] SPL {} FANOUT BEFORE: {}\n", node, last_spl, spl, static_cast<int>(ntk._storage->nodes[spl].data[0].h1));
          DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(spl));
          ntk._storage->nodes[spl].data[0].h1++;  // Increment fan-out of the splitter.
          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] SPL {} FANOUT  AFTER: {}\n", node, last_spl, spl, static_cast<int>(ntk._storage->nodes[spl].data[0].h1));
          DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(spl));

          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] NODE {} FANOUT BEFORE: {}\n", node, last_spl, node, static_cast<int>(ntk._storage->nodes[node].data[0].h1));
          DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(node));
          ntk._storage->nodes[node].data[0].h1--; // Decrement fan-out of the current node.
          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] NODE {} FANOUT  AFTER: {}\n", node, last_spl, node, static_cast<int>(ntk._storage->nodes[node].data[0].h1));
          DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(node));

          // // Notify listeners of the modification.
          // for (auto const& fn : ntk._events->on_modified)
          // {
          //   (*fn)(*it, old_preds);
          // }
          break;
        }
      }
      last_spl = spl;
      splitters.push_back(spl);
    }

    // Process the last fanout.
    // auto& preds = ntk._storage->nodes[fanouts.back()].children;
    DEBUG_PRINT("\t\t[NODE {}, LAST_SPL {}] UPDATING CONNECTIONS TO LAST FANOUT {}\n", node, last_spl, fanouts.back());
    // for (auto& pred : preds)
    for (auto pred_it = ntk._storage->nodes[fanouts.back()].children.begin();
              pred_it < ntk._storage->nodes[fanouts.back()].children.end(); pred_it++ )
    {
      if ( ntk.is_dangling( pred_it->data ) )
      {
        continue;
      }

      if (pred_it->data == node)
      {
        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] RECORDING OLD PREDS\n", node, last_spl);
        // Store the previous connections.
        std::vector<klut::signal> old_preds(ntk._storage->nodes[fanouts.back()].children.size());
        std::transform(ntk._storage->nodes[fanouts.back()].children.begin(), ntk._storage->nodes[fanouts.back()].children.end(), old_preds.begin(), [](auto c) { return c.index; });
        printVector(old_preds, 4);

        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS BEFORE\n", node, last_spl);
        for (const auto& entry : ntk._storage->nodes[fanouts.back()].children)  { DEBUG_PRINT("\t\t\t\t{}\n", entry.data); }
        pred_it->data = last_spl;                            // Replace the connection with the last splitter.
        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS AFTER\n", node, last_spl);
        for (const auto& entry : ntk._storage->nodes[fanouts.back()].children)  { DEBUG_PRINT("\t\t\t\t{}\n", entry.data); }     


        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] SPL {} FANOUT BEFORE: {}\n", node, last_spl, last_spl, static_cast<int>(ntk._storage->nodes[last_spl].data[0].h1));
        DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(last_spl));
        ntk._storage->nodes[last_spl].data[0].h1++; // Increment fan-out of the last splitter.
        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] SPL {} FANOUT  AFTER: {}\n", node, last_spl, last_spl, static_cast<int>(ntk._storage->nodes[last_spl].data[0].h1));
        DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(last_spl));

        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] NODE {} FANOUT BEFORE: {}\n", node, last_spl, node, static_cast<int>(ntk._storage->nodes[node].data[0].h1));
        DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(node));
        ntk._storage->nodes[node].data[0].h1--;     // Decrement fan-out of the current node.
        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] NODE {} FANOUT  AFTER: {}\n", node, last_spl, node, static_cast<int>(ntk._storage->nodes[node].data[0].h1));
        DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(node));

        // Notify listeners of the modification.
        for (auto const& fn : ntk._events->on_modified)
        {
          (*fn)(fanouts.back(), old_preds);
        }

        break;
      }
    }

    auto fo_ntk { mockturtle::fanout_view( ntk ) };
    for (const auto & node : fanouts)
    {
      DEBUG_PRINT("\t\t\t node: {}\n", node);
      
      fo_ntk.foreach_fanin( node, [&](const klut::signal & fi_node)
      {
        if ( fo_ntk.is_dangling( fi_node ) )
        {
          return;
        }

        DEBUG_PRINT("\t\t\t\t fanin : {}\n", fi_node);
      });
      fo_ntk.foreach_fanout( node, [&](const klut::signal & fo_node)
      {
        if ( fo_ntk.is_dangling( fo_node ) )
        {
          return;
        }

        DEBUG_PRINT("\t\t\t\t fanout: {}\n", fo_node);
      });
    }
    for (const auto & node : splitters)
    {
      DEBUG_PRINT("\t\t\t spl : {}\n", node);
      
      fo_ntk.foreach_fanin( node, [&](const klut::signal & fi_node)
      {
        if ( fo_ntk.is_dangling( fi_node ) )
        {
          return;
        }

        DEBUG_PRINT("\t\t\t\t fanin : {}\n", fi_node);
      });
      fo_ntk.foreach_fanout( node, [&](const klut::signal & fo_node)
      {
        if ( fo_ntk.is_dangling( fo_node ) )
        {
          return;
        }

        DEBUG_PRINT("\t\t\t\t fanout: {}\n", fo_node);
      });
    }
    // Ensure that the current node's fan-out count is now 1 (since all other fanouts have been replaced by splitters).
    assert( ntk._storage->nodes[node].data[0].h1 == 1 );
  } );
}

/// @brief Structure representing the potential DFF location uniquely defined by fanin, fanout and stage
struct DFF_var 
{
  klut::signal fanin;
  klut::signal fanout;
  uint32_t sigma;
  std::unordered_set<uint64_t> parent_hashes;

  DFF_var(klut::signal _fanin, klut::signal _fanout, uint32_t _sigma, std::unordered_set<uint64_t> _parent_hashes = {})
      : fanin(_fanin), fanout(_fanout), sigma(_sigma), parent_hashes(_parent_hashes) {}

  DFF_var(const DFF_var& other)
      : fanin(other.fanin), fanout(other.fanout), sigma(other.sigma), parent_hashes(other.parent_hashes) {}

  DFF_var( klut::signal _index, uint32_t _sigma, std::unordered_set<uint64_t> _parent_hashes = {} )
      : fanin( 0u ), fanout( _index ), sigma( _sigma ), parent_hashes( _parent_hashes ) {}

  DFF_var( uint64_t dff_hash, std::unordered_set<uint64_t> _parent_hashes = {} )
      : fanin( (uint64_t)dff_hash >> 40 ), fanout( (uint64_t)( dff_hash << 24 ) >> 40 ), sigma( (uint64_t)( dff_hash << 48 ) >> 40 ), parent_hashes( _parent_hashes ) {}

  std::string str() const
  {
    if ( fanin == 0 )
    {
      return fmt::format( "gate_{}_{}", fanout, sigma );
    }

    return fmt::format("var_{}_{}_{}", fanin, fanout, sigma);
  }
};

uint64_t dff_hash(klut::signal _fanin, klut::signal _fanout, uint32_t _sigma)
{
  return ( (uint64_t)_fanin << 40 ) | ( (uint64_t)_fanout << 16 ) | _sigma;
}

uint64_t dff_hash( DFF_var const& dff )
{
  return ( (uint64_t)dff.fanin << 40 ) | ( (uint64_t)dff.fanout << 16 ) | dff.sigma;
}

struct DFF_registry
{
  phmap::flat_hash_map<uint64_t, DFF_var> variables;

  DFF_var & at(node_t _fanin, node_t _fanout, glob_phase_t _sigma)
  {
    return variables.at( dff_hash(_fanin, _fanout, _sigma) );
  } 
  DFF_var & at(uint64_t _hash)
  {
    return variables.at( _hash );
  } 
  uint64_t add(node_t _fanin, node_t _fanout, glob_phase_t _phase, std::unordered_set<uint64_t> _parent_hashes = {})
  {
    uint64_t _hash = dff_hash(_phase, _fanout, _fanin);
    DFF_var temp { _fanin, _fanout, _phase, _parent_hashes };
    variables.emplace(_hash, temp);
    return _hash;
  } 

  std::string str(uint64_t hash, bool negated = false) const
  {
    const DFF_var & dff = variables.at(hash);
    if (negated)
    {
      return fmt::format( "var_{}_{}_{}.Not()", dff.fanin, dff.fanout, dff.sigma );
    }
    return fmt::format( "var_{}_{}_{}", dff.fanin, dff.fanout, dff.sigma );
  }
};


union DFF_union
{
  struct 
  {
    unsigned fanin : 24;
    unsigned fanout: 24;
    unsigned sigma : 16;
  };
  uint64_t hash;
  DFF_union(const klut::signal _fanin, const klut::signal _fanout, const uint32_t _sigma)
    : fanin(_fanin), fanout(_fanout), sigma(_sigma) {}

  DFF_union(const uint64_t _hash)
    : hash(_hash) {};

  DFF_union(const DFF_var& other)
      : fanin(other.fanin), fanout(other.fanout), sigma(other.sigma) {}

  std::string str() const
  {
      return fmt::format("var_{}_{}_{}", fanin, fanout, sigma);
  }

  // Define a custom less-than operator for DFF_union
  bool operator<(const DFF_union& other) const
  {
    if (sigma != other.sigma) return sigma < other.sigma;
    if (fanin != other.fanin) return fanin < other.fanin;
    return fanout < other.fanout;
  }
};

// Define a custom hash function for DFF_union
struct DFF_unionHash
{
  std::size_t operator()(const DFF_union& dff) const
  {
    return static_cast<std::size_t>(dff.hash);
  }
};

std::tuple<
  std::vector<DFF_union>,            // helper_vars: always true 
  std::vector<DFF_union>,            // sa_dff : always true + contribute to cost
  std::set<std::vector<DFF_union>>,  // stage_constraints: no more than one DFF per sigma
  std::set<std::deque<DFF_union>>,   // buffer_constraints: at least one buffer per n_phases
  std::set<std::vector<DFF_union>>,  // inverted_t1_input: if AS/SA/T1 -> T1 and inverted, need 1+ DFF 
  std::set<std::vector<DFF_union>>,   // merger_t1_input : if merger -> T1, need 1+ DFF (avoid double pulse)
  std::map<klut::signal, std::set<std::vector<DFF_union>>> //truncated_t1_paths: extract hashes and 
> dff_from_threads(const klut & ntk, const Path & path, const uint8_t & n_phases, const phmap::flat_hash_map<klut::node, std::array<bool, 3>> & input_phases, bool verbose = true)
{
  // always true vars
  std::vector<DFF_union> helper_vars;
  // SA DFF
  std::vector<DFF_union> sa_dff;
  // no more than one DFF per sigma
  std::set<std::vector<DFF_union>>  stage_constraints; 
  // at least one buffer per n_phases
  std::set<std::deque<DFF_union>>                 buffer_constraints;
  // if AS/SA/T1 -> T1 and inverted, at least one DFF should be placed (will be converted to an inverter)
  // no need to create constraint if the distance is greater than n_phases
  std::set<std::vector<DFF_union>>   inverted_t1_input;
  // if merger -> T1, at least one DFF should be placed (avoid double pulse)
  // no need to create constraint if the distance is greater than n_phases
  std::set<std::vector<DFF_union>>     merger_t1_input;
  // the last DFF for each T1 input should have different phase
  // create constraint until n_phases from the T1 cell
  // std::vector<std::vector<std::vector<DFF_union>>>   unequal_t1_inputs;

  // for each T1 gate (klut::signal), extract three truncated paths (set)
  // each path is the vector of DFF_var hashes <std::vector<DFF_union>>
  std::map<klut::signal, std::set<std::vector<DFF_union>>> truncated_t1_paths;

  const std::vector<std::vector<klut::signal>> path_threads = path.path_threads(ntk, verbose);
  std::vector<std::map<uint32_t, std::vector<DFF_union>>> thread_dffs_by_sigma;
  
  for (const std::vector<klut::signal> & path_thread : path_threads)
  {
    const klut::signal & path_target = path_thread.front();
    NodeData target_data { ntk.value( path_target ) }; 

    DEBUG_PRINT("[DFF F T] PROCESSING {{ {} }}\n", fmt::join(path_thread, "<-"));
    DEBUG_PRINT("[DFF F T] FRONT: {} TYPE: {} SIGMA: {}\n", path_target, GATE_TYPE[target_data.type], static_cast<int>(target_data.sigma));

    // First, generate DFF variables
    std::map<uint32_t, std::vector<DFF_union>> dff_by_sigma;
    // iterating until the source
    for (auto it = path_thread.begin(); it != path_thread.end()-1; ++it)
    {
      const klut::signal & fo_node = *it;
      const NodeData fo_data { ntk.value(fo_node) };
      const klut::signal & fi_node = *(it + 1);
      const NodeData fi_data { ntk.value(fi_node) };
      
      DEBUG_PRINT("\t between {} ({},{}) and {} ({},{}) \n", fo_node, GATE_TYPE[fo_data.type], static_cast<int>(fo_data.sigma), fi_node, GATE_TYPE[fi_data.type], static_cast<int>(fi_data.sigma));
      
      uint32_t max_sigma = fo_data.sigma - static_cast<int>(fo_data.type == AS_GATE || fo_data.type == T1_GATE);
      if (max_sigma == target_data.sigma) 
      { 
        max_sigma--; 
      }
      for (auto sigma = fi_data.sigma; sigma <= max_sigma; ++sigma)
      {
        dff_by_sigma[sigma].emplace_back( fi_node, fo_node, sigma );
        DEBUG_PRINT("\t Adding {} \n", dff_by_sigma[sigma].back().str());
      }

      if (fo_data.type == SA_GATE)
      {
        // do not create if the preceding gate is AS and at the same sigma
        if (!(fi_data.type == AS_GATE && fi_data.sigma == fo_data.sigma))
        {
          sa_dff.emplace_back( fi_node, fo_node, fo_data.sigma );
        }
      }
    }
    // create DFF_var for the source gate (always true)
    const klut::signal & fo_node = path_thread.back();
    NodeData fo_data { ntk.value(fo_node) };

    const DFF_union dff { 0, fo_node, fo_data.sigma };
    dff_by_sigma[fo_data.sigma].push_back(dff);
    helper_vars.push_back(dff);


    // stage & buffer constraints
    std::deque<DFF_union> chain;
    for (const auto & [sigma, var_hashes] : dff_by_sigma)
    {
      DEBUG_PRINT("[SIGMA={}] ", sigma); 
      if (verbose)
      {
        for (const DFF_union & dff : var_hashes)
        {
          DEBUG_PRINT(" {} ", dff.str()); 
        }
      }
      DEBUG_PRINT("\n"); 

      if (var_hashes.size() > 1)
      {
        stage_constraints.emplace( var_hashes );
        DEBUG_PRINT("\tCreated stage constraint\n"); 
      }
      
      // insert new variables at the back
      chain.insert(chain.end(), var_hashes.begin(), var_hashes.end());
      const int min_sigma = sigma - n_phases;
      // remove old variables from the front
      while ( sigma - chain.front().sigma >= n_phases)
      {
        chain.pop_front();
      }

      DEBUG_PRINT("[CHAIN] "); 
      if (verbose)
      {
        for (const DFF_union & dff : chain)
        {
          DEBUG_PRINT(" {} ", dff.str()); 
        }
      }
      DEBUG_PRINT("\n"); 

      DEBUG_PRINT("Chain span : {} - {} = {} (n_phases = {})\n", 
        static_cast<int>(chain.back().sigma),
        static_cast<int>(chain.front().sigma), 
        static_cast<int>(chain.back().sigma, chain.front().sigma),
        n_phases 
        ); 

      if ( chain.back().sigma - chain.front().sigma == (n_phases - 1) )
      {
        DEBUG_PRINT("\t[CHAIN] emplaced\n"); 
        buffer_constraints.emplace(chain);
      }
    }

    if ( target_data.type == T1_GATE )
    {
      DEBUG_PRINT("Target is T1\n"); 
      std::vector<klut::signal> fanins;

      // check input complementation
      ntk.foreach_valid_fanin( path_target, [&]( const klut::signal & fanin ) { fanins.push_back( fanin ); });
      uint8_t fanin_ID = std::find(fanins.begin(), fanins.end(), path_thread[1]) - fanins.begin();
      DEBUG_PRINT("\t Fanin ID {}\n", fanin_ID); 
      assert(fanin_ID < 3);

      bool is_negated = input_phases.at(path_target)[fanin_ID];
      DEBUG_PRINT("\t Fanin is {}negated\n", is_negated?"":"NOT "); 

      // loop until merger or until the source
      std::vector<DFF_union> dff_chain;
      for (auto it = path_thread.begin(); it != path_thread.end()-1; ++it)
      {
        const klut::signal & fo_node = *it;
        NodeData fo_data { ntk.value(fo_node) };
        const klut::signal & fi_node = *(it + 1);
        NodeData fi_data { ntk.value(fi_node) };

        DEBUG_PRINT("\t [t1] between {} ({},{}) and {} ({},{}) \n", fo_node, GATE_TYPE[fo_data.type], static_cast<int>(fo_data.sigma), fi_node, GATE_TYPE[fi_data.type], static_cast<int>(fi_data.sigma));

        uint8_t fi_size = 0u;
        ntk.foreach_valid_fanin(fi_node, [&fi_size](const klut::signal & fi_node) { fi_size++; });
        DEBUG_PRINT("\t [t1] fanin size: {} \n", fi_size);

        assert(fo_data.type != AS_GATE && fo_data.type != SA_GATE);

        auto max_sigma = fo_data.sigma - static_cast<int>( fo_data.type == T1_GATE );
        max_sigma = std::min(max_sigma, target_data.sigma - 1);
        for (auto sigma = fi_data.sigma; sigma <= max_sigma; ++sigma)
        {
          DEBUG_PRINT("\t [t1] accessing {} ({},{}) \n", fi_node, fo_node, sigma);
          dff_chain.emplace_back( fi_node, fo_node, sigma );
          DEBUG_PRINT("\t [t1] pushing {} \n", dff_chain.back().str());
        }
        if ( fi_data.type == AA_GATE )
        {
          // found splitter, check the next gate
          if (fi_size == 1) 
          { 
            DEBUG_PRINT("\t [t1] found splitter \n");
            continue; 
          }
          // found a merger, record the chain
          else if (fi_size == 2)
          {
            DEBUG_PRINT("\t [t1] found merger \n");

            // check whether the chain is short enough
            if ( target_data.sigma - fi_data.sigma <= n_phases )
            {
              DEBUG_PRINT("\t [t1] Recording the merger chain ");
              merger_t1_input.emplace( dff_chain );
            }

            DEBUG_PRINT("\t [t1] Full t1 chain: ");
            for (const DFF_union &  dff : dff_chain)
            {
              DEBUG_PRINT(" {}", dff.str());
            }
            DEBUG_PRINT("\n");

            DEBUG_PRINT("\t [t1] Recording the t1 chain ");
            const auto min_sigma = target_data.sigma - n_phases;
            DEBUG_PRINT("\t [t1] min_sigma =  {} - {} = {}\n", static_cast<int>(target_data.sigma), n_phases, min_sigma);
            std::vector<DFF_union> truncated_chain;
            for (const DFF_union &  dff : dff_chain)
            {
              if (dff.sigma >= min_sigma)
              {
                truncated_chain.push_back(dff);
                DEBUG_PRINT(" {}", dff.str());
              }
            }
            truncated_t1_paths[path_target].emplace( truncated_chain );
            DEBUG_PRINT("\n");
            break;
          }
          else //2+ input merger ???
          {
            throw;
          }
        }
        // adder, bar, max, 

        // no AA -> T1 connection
        // rather, AS/SA/T1 -> T1 connection
        // add to the truncated_t1_paths (can add the entire path)
        else
        {
          // Finish the chain with the gate variable
          DEBUG_PRINT("\t [t1] accessing {} ({},{}) \n", 0, fi_node, static_cast<int>(fi_data.sigma));
          dff_chain.emplace_back( 0, fi_node, static_cast<int>(fi_data.sigma) );
          DEBUG_PRINT("\t [t1] finishing the chain {} \n", dff_chain.back().str());
        
          if (is_negated)
          {
            DEBUG_PRINT("\t [t1] negated AS/SA/T1 -> T1 connection \n");
            
            // check whether the chain is short enough
            if ( target_data.sigma - fi_data.sigma <= n_phases )
            {
              DEBUG_PRINT("\t [t1] Recording the negated chain ");
              inverted_t1_input.emplace( dff_chain );
            }
          }
          else 
          {
            DEBUG_PRINT("\t [t1] positive AS/SA/T1 -> T1 connection \n");
          }

          DEBUG_PRINT("\t [t1] Full t1 chain: ");
          for (const DFF_union &  dff : dff_chain)
          {
            DEBUG_PRINT(" {}", dff.str());
          }
          DEBUG_PRINT("\n");

          DEBUG_PRINT("\t [t1] Recording the t1 chain ");
          const auto min_sigma = target_data.sigma - n_phases;
          DEBUG_PRINT("\t [t1] min_sigma =  {} - {} = {}\n", static_cast<int>(target_data.sigma), n_phases, min_sigma);
          std::vector<DFF_union> truncated_chain;
          for (const DFF_union &  dff : dff_chain)
          {
            if (dff.sigma >= min_sigma)
            {
              truncated_chain.push_back(dff);
              DEBUG_PRINT(" {}", dff.str());
            }
          }
          truncated_t1_paths[path_target].emplace( truncated_chain );
          DEBUG_PRINT("\n");
          break;
        }
      }
    }
  }
  return std::make_tuple( helper_vars, sa_dff, stage_constraints, buffer_constraints, inverted_t1_input, merger_t1_input, truncated_t1_paths );
}

void write_dff_cfg(
  const klut & ntk,
  const std::string & filename,
  const std::vector<DFF_union> & helper_vars,            // helper_vars: always true 
  const std::vector<DFF_union> & sa_dff,            // helper_vars: always true 
  const std::set<std::vector<DFF_union>> & stage_constraints,  // stage_constraints: no more than one DFF per sigma
  const std::set<std::deque<DFF_union>> & buffer_constraints,   // buffer_constraints: at least one buffer per n_phases
  const std::set<std::vector<DFF_union>> & inverted_t1_input,  // inverted_t1_input: if AS/SA/T1 -> T1 and inverted, need 1+ DFF 
  const std::set<std::vector<DFF_union>> & merger_t1_input,   // merger_t1_input : if merger -> T1, need 1+ DFF (avoid double pulse)
  const std::map<klut::signal, std::set<std::vector<DFF_union>>> & truncated_t1_paths,
  bool verbose = false
  )
{
  // auto out = fmt::output_file(filename);
  std::ofstream os { filename };
  // always true variables, do not contribute to cost
  for (const DFF_union dff : helper_vars)
  {
    DEBUG_PRINT("\t [helper_vars] accessing dff: {} \n", dff.str());
    os << fmt::format("HELPER,{}\n", dff.str());
  }
  // always true variables,  contribute to cost
  for (const DFF_union dff : sa_dff)
  {
    DEBUG_PRINT("\t [sa_dff] accessing dff: {} \n", dff.str());
    os << fmt::format("SA_REQUIRED,{}\n", dff.str());
  }
  // stage_constraints: no more than one DFF per sigma
  for (const std::vector<DFF_union> & hashes : stage_constraints)
  {
    os << "PHASE";
    for (const DFF_union dff : hashes)
    {
      DEBUG_PRINT("\t [stage_constraints] accessing dff: {} \n", dff.str());
      os << fmt::format(",{}", dff.str());
    }
    os << fmt::format("\n");
  }
  // buffer_constraints: at least one buffer per n_phases
  for (const std::deque<DFF_union> & hashes : buffer_constraints)
  {
    os << "BUFFER";
    for (const DFF_union dff : hashes) 
    {
      DEBUG_PRINT("\t [buffer_constraints] accessing dff: {} \n", dff.str());
      os << fmt::format(",{}", dff.str());
    }
    os << fmt::format("\n");
  }
  // inverted_t1_input: if AS/SA/T1 -> T1 and inverted, need 1+ DFF 
  for (const std::vector<DFF_union> & hashes : inverted_t1_input)
  {
    os << "INVERTED_INPUT";
    for (const DFF_union dff : hashes) 
    {
      DEBUG_PRINT("\t [inverted_t1_input] accessing dff: {} \n", dff.str());
      os << fmt::format(",{}", dff.str());
    }
    os << fmt::format("\n");
  }
  // merger_t1_input: merger -> T1, need 1+ DFF
  for (const std::vector<DFF_union> & hashes : merger_t1_input)
  {
    os << "AT_LEAST_ONE";
    for (const DFF_union dff : hashes) 
    {
      DEBUG_PRINT("\t [merger_t1_input] accessing dff: {} \n", dff.str());
      os << fmt::format(",{}", dff.str());
    }
    os << fmt::format("\n");
  }
  // const std::map<klut::signal, std::set<std::vector<uint64_t>>> & truncated_t1_paths
  // merger_t1_input: merger -> T1, need 1+ DFF
  for (const auto & [t1_node, paths] : truncated_t1_paths)
  {
    assert(NodeData(ntk.value(t1_node)).type == T1_GATE);
    assert(paths.size() == 3);
    for (const std::vector<DFF_union> & path : paths)
    {
      os << fmt::format("T1_FANIN,{}",t1_node);
      for (const DFF_union & dff : path)
      {
        DEBUG_PRINT("\t [truncated_t1_paths] accessing dff: {} \n", dff.str());
        os << fmt::format(",{}", dff.str());
      }
      os << fmt::format("\n");
    }
  }
}

std::vector<klut::signal> get_fanins(const klut & ntk, const klut::signal & fo_node, const phmap::flat_hash_map<klut::signal, std::vector<klut::signal>> & cut_leaves)
{
  if (NodeData(ntk.value(fo_node)).type == T1_GATE)
  {
    // get leaves of the 3-cut
    return cut_leaves.at(fo_node);
  }
  else
  {
    std::vector<klut::signal> fanins;
    ntk.foreach_fanin(fo_node, [&](const klut::signal & fi_node)
    {
      fanins.push_back(fi_node);
    });
  }
}

/* the version of 'get_fanins' that does not distinguish T1 gates from others */
std::vector<klut::signal> get_fanins( klut const& ntk, klut::signal const& n )
{
  std::vector<klut::signal> fanins;
  ntk.foreach_fanin( n, [&]( klut::signal const& ni )
  {
    fanins.push_back( ni );
  } );
}

void write_klut_specs(const klut & ntk, const std::string & filename)
{
  auto ntk_fo = mockturtle::fanout_view<klut, false>( ntk );
  std::ofstream spec_file (filename);
  ntk.foreach_node( [&]( const klut::signal & sig ) 
  {
    std::vector<klut::node> fanouts;
    ntk_fo.foreach_fanout(sig, [&](const klut::node fo) { fanouts.push_back(fo);  });
    std::vector<klut::node> fanins;
    ntk_fo.foreach_fanin (sig, [&](const klut::node fi) {  fanins.push_back(fi);  });

    // TODO: Infer from the ntk, not Primitive object
    const NodeData & node_data {ntk.value( sig )};
    uint8_t prim_type = ( ntk.fanin_size( sig ) == 0 ) ? PI_GATE : node_data.type;
    spec_file << fmt::format("{0},{1},{2},{3}\n", sig, prim_type, fmt::join(fanins, "|"), fmt::join(fanouts, "|"));
  });
}

void write_klut_specs_T1(const klut & ntk, const phmap::flat_hash_map<klut::signal, std::vector<klut::signal>> & cut_leaves, const phmap::flat_hash_map<klut::signal, klut::signal> & representatives, const std::string & filename)
{
  auto ntk_fo = mockturtle::fanout_view<klut, false>( ntk );
  std::ofstream spec_file (filename);
  ntk.foreach_node( [&]( const klut::signal & sig ) 
  {
    if (NodeData(ntk.value(sig)).type == T1_GATE)
    {
      if (representatives.at(sig) != sig)
      {
        return;
      }
      // TODO: also need to record the bound nodes
      // need to ensure the sigmas given to T1 outputs are equal
    }
    std::vector<klut::node> fanouts;
    ntk_fo.foreach_fanout(sig, [&](const klut::node fo) { fanouts.push_back(fo);  });
    std::vector<klut::node> fanins = get_fanins(ntk, sig, cut_leaves);
    // ntk_fo.foreach_fanin (sig, [&](const klut::node fi) {  fanins.push_back(fi);  });

    const NodeData & node_data {ntk.value( sig )};
    uint8_t prim_type = ( ntk.fanin_size( sig ) == 0 ) ? PI_GATE : node_data.type;
    spec_file << fmt::format("{0},{1},{2}\n", sig, prim_type, fmt::join(fanins, "|"));
  });
}


std::vector<Path> extract_paths(const klut & ntk, bool verbose = false)
{
  DEBUG_PRINT("\t[i] ENTERED FUNCTION extract_paths\n");
  std::vector<Path> paths;

  ntk.foreach_node([&](const klut::signal & fo_node)
  {
    DEBUG_PRINT("\t\t[i] PROCESSING NODE {}\n", fo_node);
    if (ntk.is_constant(fo_node) || ntk.is_pi(fo_node))
    {
      DEBUG_PRINT("\t\t\t[NODE {}] the node is CONSTANT/PI\n", fo_node);
      return;
    }
    NodeData fo_node_data = ntk.value(fo_node);
    if ( fo_node_data.type == AA_GATE ) 
    {
      DEBUG_PRINT("\t\t\t[NODE {}] the node is AA, skipping\n", fo_node);
      return;
    }

    // at this point, the node should be AS/SA/T1
    DEBUG_PRINT("\t\t[NODE {}] the node is AS/SA/T1, continuing...\n", fo_node);

    ntk.foreach_fanin(fo_node, [&](const klut::signal & fi_node)
    {
      DEBUG_PRINT("\t\t\t[NODE {}] processing fanin {}\n", fo_node, fi_node);
      // Create a path object with only a target
      Path node_path;
      node_path.targets.emplace( fo_node );

      std::vector<klut::signal> stack { fi_node };
      
      DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}] created stack\n", fo_node, fi_node);
      
      std::set<klut::signal> seen;
      while (!stack.empty())
      {
        DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}] stack contents:\n", fo_node, fi_node);
        printVector(stack, 4);

        const klut::signal & n = stack.back();
        stack.pop_back();

        DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}]: Analyzing node {}\n", fo_node, fi_node, n);

        // A constant does not have any effect on the DFF placement, we can skip it
        if ( ntk.is_constant( n ) )
        {
          continue;
        }        
        const NodeData n_data { ntk.value(n) };
        // Found a source of the path, add to sources, do not continue traversal
        if ( ntk.is_pi(n) || n_data.type == AS_GATE || n_data.type == SA_GATE || n_data.type == T1_GATE )
        {
          DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}]: node {} is a source \n", fo_node, fi_node, n);
          node_path.sources.emplace( n );
        }
        // Found AA gate, add to internal nodes, add parents for further traversal
        else if ( n_data.type == AA_GATE )
        {
          DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}]: node is INTERNAL adding fanins \n", fo_node, fi_node, n);
          node_path.internals.emplace( n );

          ntk.foreach_fanin(n, [&](const klut::signal & sig){
            stack.push_back( sig );
          });
        }
        else
        {
          DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}]: Signal {}: {} is not recognized \n", fo_node, fi_node, n, GATE_TYPE.at( n_data.type ));
          throw "Unsupported case";
        }
        seen.emplace( n );
      }

      // Identify overlapping paths
      std::vector<size_t> to_merge;
      for (size_t i = 0u; i < paths.size(); ++i)
      {
        Path & known_paths = paths[i];
        // merge if there are sources in common
        if( haveCommonElements( known_paths.sources, node_path.sources) )
        {
          to_merge.push_back(i);
        }
      }
      // Merge overlapping paths into the path object and remove the path object
      // Iterating in reverse order to preserve order
      for (auto it = to_merge.rbegin(); it != to_merge.rend(); ++it) 
      {
        auto idx = *it;
        // fmt::print("Before absorption\n");
        // node_path.print();
        node_path.absorb(paths[idx]);
        // fmt::print("After absorption\n");
        // node_path.print();
        paths.erase(paths.begin() + idx);
      }
      paths.push_back(node_path);
    });
  });
  return paths;
}

void merge_overlapping_paths( std::vector<Path>& paths, Path& current_path )
{
  /* identify overlapping paths that can be merged */
  std::vector<size_t> index_to_merge;
  for ( auto i{ 0u }; i < paths.size(); ++i )
  {
    Path& known_path = paths[i];
    if ( haveCommonElements( known_path.sources, current_path.sources ) )
    {
      index_to_merge.push_back( i );
    }
  }

  /* merge identified overlapping paths */
  for ( auto it{ index_to_merge.rbegin() }; it != index_to_merge.rend(); ++it )
  {
    auto idx = *it;
    current_path.absorb( paths[idx] );
    paths.erase( paths.begin() + idx );
  }
}

void buildPath(const klut & ntk, Path & node_path, std::vector<klut::signal> stack, const bool verbose = false)
{
  std::set<klut::signal> seen;
  while (!stack.empty())
  {
    DEBUG_PRINT("\t\t\tStack contents:\n");
    printVector(stack, 4);

    const klut::signal & n = stack.back();
    stack.pop_back();

    DEBUG_PRINT("\t\t\tAnalyzing node {}\n", n);

    // A constant does not have any effect on the DFF placement, we can skip it
    if ( ntk.is_constant( n ) )
    {
      continue;
    }        
    const NodeData n_data { ntk.value(n) };
    // Found a source of the path, add to sources, do not continue traversal
    if ( ntk.is_pi(n) || n_data.type == AS_GATE || n_data.type == SA_GATE || n_data.type == T1_GATE )
    {
      DEBUG_PRINT("\t\t\tnode {} is a source \n", n);
      node_path.sources.emplace( n );
    }
    // Found AA gate, add to internal nodes, add parents for further traversal
    else if ( n_data.type == AA_GATE )
    {
      DEBUG_PRINT("\t\t\tnode is INTERNAL adding fanins \n", n);
      node_path.internals.emplace( n );

      ntk.foreach_fanin(n, [&](const klut::signal & sig){
        stack.push_back( sig );
      });
    }
    else
    {
      DEBUG_PRINT("\t\t\tSignal {}: {} is not recognized \n", n, GATE_TYPE.at( n_data.type ));
      throw "Unsupported case";
    }
    seen.emplace( n );
  }
}

/// @brief same as extract_paths but with the support of T1 cells. Uses helper functions
/// @param ntk 
/// @param verbose 
/// @return 
std::vector<Path> extract_paths_t1(const klut & ntk, const phmap::flat_hash_map<klut::signal, klut::signal> & representatives, bool verbose = false)
{
  DEBUG_PRINT("\t[i] ENTERED FUNCTION extract_paths\n");
  std::vector<Path> paths;

  mockturtle::topo_view<klut>( ntk ).foreach_node( [&](const klut::signal & fo_node)
  {
    DEBUG_PRINT("\t\t[i] PROCESSING NODE {}\n", fo_node);
    if (ntk.is_constant(fo_node) || ntk.is_pi(fo_node))
    {
      DEBUG_PRINT("\t\t\t[NODE {}] the node is PI/CONSTANT\n", fo_node);
      return;
    }
    if ( ntk.is_dangling( fo_node ) ) 
    {
      DEBUG_PRINT("\t\t\t[NODE {}] the node is dangling, skipping\n", fo_node);
      return;
    }
    NodeData fo_node_data = ntk.value(fo_node);
    if ( fo_node_data.type == AA_GATE ) 
    {
      DEBUG_PRINT("\t\t\t[NODE {}] the node is AA, skipping\n", fo_node);
      return;
    }

    // at this point, the node should be AS/SA/T1
    DEBUG_PRINT("\t\t[NODE {}] the node is AS/SA/T1, continuing...\n", fo_node);

    // If the gate is T1:
    //   - all leaves of the cut belong to the same path due to the transverse constraint (unlike AS/SA)
    //   - only a single node among T1 outputs needs to be considered. Other outputs would represent the same path
    if (fo_node_data.type == T1_GATE)
    {
      // Check if the node is the main representative of the T1 cell
      if (representatives.at(fo_node) != fo_node)
      {
        // if not, skip, the T1 cell was (or will be) considered by another node
        return;
      }

      Path node_path;
      node_path.targets.emplace( fo_node );
      std::vector<klut::signal> stack;
      ntk.foreach_valid_fanin( fo_node, [&]( auto const& ni ) {
        stack.push_back( ni );
      } );
      buildPath(ntk, node_path, stack, verbose);
      merge_overlapping_paths( paths, node_path );
      paths.push_back(node_path);
    }
    else if (fo_node_data.type == AS_GATE || fo_node_data.type == SA_GATE)
    {
      Path node_path;
      node_path.targets.emplace( fo_node );
      ntk.foreach_valid_fanin( fo_node, [&] (const klut::signal & fi_node) {
        std::vector<klut::signal> stack { fi_node };
        buildPath(ntk, node_path, stack, verbose);
        merge_overlapping_paths( paths, node_path );
        paths.push_back(node_path);
      } );
    }
    else
    {
      DEBUG_PRINT("\t\t\tSignal {}: {} is not recognized \n", fo_node, GATE_TYPE.at( fo_node_data.type ));
      throw "Unsupported case";
    }
  } );

  return paths;
}

/// @brief Create binary variables for DFF placement in a given path
/// @param path - a path object to insert DFFs into
/// @param NR - unordered_map of NtkNode objects 
/// @param n_phases - # of phases
/// @param verbose - prints debug messages if set to *true*
/// @return 
std::tuple<DFF_registry, uint64_t, std::vector<uint64_t>> dff_vars_single_paths(const Path & path, const klut & ntk, const uint8_t n_phases, bool verbose = false)
{
  DFF_registry DFF_REG;
  std::vector<uint64_t> required_SA_DFFs;

  std::vector<std::tuple<klut::signal, uint64_t>> stack;
  for (const klut::signal & tgt : path.targets)
  {
    stack.emplace_back(tgt, 0);
  }
  DEBUG_PRINT("[DFF] Target nodes: {}\n", fmt::join(path.targets, ","));

  auto precalc_ndff = 0u;
  
  while (!stack.empty())
  { 
    if (verbose)
    {
      DEBUG_PRINT("STACK :\n");
      for (const auto & [ fo_node, earliest_child_hash ] : stack)
      {
        NodeData fo_data { ntk.value( fo_node ) };
        if (earliest_child_hash != 0)
        {
          DEBUG_PRINT("\t{}({})[{}], {}\n", GATE_TYPE.at((int)fo_data.type), fo_node, (int)fo_data.sigma, DFF_REG.at(earliest_child_hash).str());
        }
        else
        {
          DEBUG_PRINT("\t{}({})[{}]\n", GATE_TYPE.at((int)fo_data.type), fo_node, (int)fo_data.sigma);
        }
      }
    }
    // fixing the C++17 bug with structured bindings
    auto _node_tuple = stack.back();
    const klut::signal fo_node          = std::get<0>(_node_tuple);
    const uint64_t earliest_child_hash  = std::get<1>(_node_tuple);
    stack.pop_back();
    NodeData fo_data { ntk.value( fo_node ) };

    // AS gates and T1 gate are clocked, so one needs to start one stage earlier
    uint32_t latest_sigma = fo_data.sigma - (fo_data.type == AS_GATE || fo_data.type == T1_GATE);
    DEBUG_PRINT("[DFF] Analyzing child: {}({})[{}]\n", GATE_TYPE.at(fo_data.type), fo_node, (int)fo_data.sigma);

    ntk.foreach_fanin(fo_node, [&](const klut::signal & fi_node)
    {
      NodeData fi_data { ntk.value( fi_node ) };
      /* if fi_node is an AA gate, it is allowed to put a DFF at the phase that fi_node is assigned to */
      uint32_t earliest_sigma = fi_data.sigma + (fi_data.type != AA_GATE);

      DEBUG_PRINT("\t[DFF] Analyzing parent: {}({})[{}]\n", GATE_TYPE.at(fi_data.type), fi_node, (int)fi_data.sigma);

      // check if the chain is straight - #DFF is just floor(delta-phase), no need to create the dff vars
      if (fo_data.type != AA_GATE && fi_data.type != AA_GATE)
      {
        // special case when an AS gate feeds directly into SA gate
        if (fo_data.sigma == fi_data.sigma)
        {
          DEBUG_PRINT("\t[DFF] Straight chain: AS{} -> SA{}\n", fi_node, fo_node);
          // do nothing, no additional DFFs needed
          assert(fo_data.type == SA_GATE && fi_data.type == AS_GATE && ntk.fanout_size(fi_node) == 1);
        }
        else
        {
          DEBUG_PRINT("\t[DFF] Straight chain: {}[{}] -> {}[{}]\n", GATE_TYPE.at(fi_data.type), (int)fi_data.sigma, GATE_TYPE.at(fo_data.type), (int)fo_data.sigma);
          // straight chain, just floor the difference!
          precalc_ndff += (fo_data.sigma - fi_data.sigma - 1)/n_phases + (fo_data.type == SA_GATE); //extra DFF before SA gate
        }
        return;
      }

      DEBUG_PRINT("\t[DFF] Non-straight chain: {}[{}] -> {}[{}]\n", GATE_TYPE.at(fi_data.type), (int)fi_data.sigma, GATE_TYPE.at(fo_data.type), (int)fo_data.sigma);
      std::vector<uint64_t> out_hashes;
      DEBUG_PRINT("\tAdding new DFFs [reg size = {}]\n", DFF_REG.variables.size());

      for (glob_phase_t sigma = earliest_sigma; sigma <= latest_sigma; ++sigma)
      {
        uint64_t new_hash = DFF_REG.add(fi_node, fo_node, sigma);
        out_hashes.push_back(new_hash);
        DEBUG_PRINT("\tAdded new DFFs at phase {} [reg size = {}]\n", sigma, DFF_REG.variables.size());
      }
      DEBUG_PRINT("\tConnecting new DFFs\n");
      for (auto i = 1u; i < out_hashes.size(); ++i)
      {
        DFF_var & dff = DFF_REG.at( out_hashes[i] );
        dff.parent_hashes.emplace(out_hashes[i-1]);
      }
      if (fo_data.type == SA_GATE)
      {
        // ensure that the SA gate is placed at least one stage after the fanin 
        //(to leave space for a DFF)
        assert( !out_hashes.empty() );
        // The last DFF in the chain is required
        required_SA_DFFs.push_back(out_hashes.back());
      }
      // if there are DFFs, the earliest_hash is the first hash in the chain
      // otherwise, it is the earliest_hash of the previous chain
      uint64_t earliest_hash = (out_hashes.empty()) ? earliest_child_hash : out_hashes.front();
      // if the node is internal, connect with the fanout phase
      if (fo_data.type == AA_GATE && !out_hashes.empty() && earliest_hash != 0 && earliest_child_hash != 0)
      {
        DFF_var & child_dff = DFF_REG.at( earliest_child_hash );
        DEBUG_PRINT("\tPrior node is {}[{}]\n", child_dff.str(), (int)child_dff.sigma); 
        // assert(child_dff.fanin == fo_id);
        child_dff.parent_hashes.emplace( out_hashes.back() );
      }
      if (fi_data.type == AA_GATE)
      {
        stack.emplace_back( fi_node, earliest_hash );
        DEBUG_PRINT("\tEmplacing {}({})[{}], {}\n", GATE_TYPE.at(fi_data.type), fi_node, (int)fi_data.sigma, (earliest_hash!=0)?DFF_REG.at(earliest_hash).str():"");
      }
    });
  }
  return std::make_tuple(DFF_REG, precalc_ndff, required_SA_DFFs);
}

/// @brief Create binary variables for DFF placement in a given path (with support of T1 cells)
/// @param path - a path object to insert DFFs into
/// @param NR - unordered_map of NtkNode objects 
/// @param n_phases - # of phases
/// @param verbose - prints debug messages if set to *true*
/// @return 
std::tuple<DFF_registry, uint64_t, std::vector<uint64_t>> dff_vars_single_paths_t1( const Path & path, const klut & ntk, const uint8_t n_phases, 
                                                                                    phmap::flat_hash_map<klut::node, std::array<uint64_t, 3>>& DFF_closest_to_t1s, 
                                                                                    bool verbose = false )
{
  DFF_registry DFF_REG;
  std::vector<uint64_t> required_SA_DFFs;

  std::vector<std::tuple<klut::signal, uint64_t>> stack;
  for (const klut::signal & tgt : path.targets)
  {
    stack.emplace_back(tgt, 0);
  }
  DEBUG_PRINT("[DFF] Target nodes: {}\n", fmt::join(path.targets, ","));

  auto precalc_ndff = 0u;
  
  while (!stack.empty())
  { 
    if (verbose)
    {
      DEBUG_PRINT("STACK :\n");
      for (const auto & [ fo_node, earliest_child_hash ] : stack)
      {
        NodeData fo_data { ntk.value( fo_node ) };
        if (earliest_child_hash != 0)
        {
          DEBUG_PRINT("\t{}({})[{}], {}\n", GATE_TYPE.at((int)fo_data.type), fo_node, (int)fo_data.sigma, DFF_REG.at(earliest_child_hash).str());
        }
        else
        {
          DEBUG_PRINT("\t{}({})[{}]\n", GATE_TYPE.at((int)fo_data.type), fo_node, (int)fo_data.sigma);
        }
      }
    }
    // fixing the C++17 bug with structured bindings
    auto _node_tuple = stack.back();
    const klut::signal fo_node          = std::get<0>(_node_tuple);
    const uint64_t earliest_child_hash  = std::get<1>(_node_tuple);
    stack.pop_back();
    NodeData fo_data { ntk.value( fo_node ) };

    // AS gates and T1 gate are clocked, so one needs to start one stage earlier
    uint32_t latest_sigma = fo_data.sigma - (fo_data.type == AS_GATE || fo_data.type == T1_GATE);
    DEBUG_PRINT("[DFF] Analyzing child: {}({})[{}]\n", GATE_TYPE.at(fo_data.type), fo_node, (int)fo_data.sigma);
    bool is_t1{ fo_data.type == T1_GATE };
    std::array<uint64_t, 3> DFF_closest_to_t1;
    uint8_t fanin_ind{ 0u };

    ntk.foreach_valid_fanin(fo_node, [&](const klut::signal & fi_node)
    {
      NodeData fi_data { ntk.value( fi_node ) };
      uint32_t earliest_sigma = fi_data.sigma + static_cast<uint32_t>( fi_data.type != AA_GATE);

      DEBUG_PRINT("\t[DFF] Analyzing parent: {}({})[{}]\n", GATE_TYPE.at(fi_data.type), fi_node, (int)fi_data.sigma);

      // check if the chain is straight - #DFF is just floor(delta-phase), no need to create the DFF vars
      if (  (fo_data.type == AS_GATE || fo_data.type == SA_GATE) // the gate should be AS or SA
            && fi_data.type != AA_GATE )                         // the fanin to the gate should be AS/SA/T1
      {
        // special case when an AS gate feeds directly into SA gate
        if (fo_data.sigma == fi_data.sigma)
        {
          DEBUG_PRINT("\t[DFF] Straight chain: AS{} -> SA{}\n", fi_node, fo_node);
          // do nothing, no additional DFFs needed
          assert(fo_data.type == SA_GATE && fi_data.type == AS_GATE && ntk.fanout_size(fi_node) == 1);
        }
        else
        {
          DEBUG_PRINT("\t[DFF] Straight chain: {}[{}] -> {}[{}]\n", GATE_TYPE.at(fi_data.type), (int)fi_data.sigma, GATE_TYPE.at(fo_data.type), (int)fo_data.sigma);
          // straight chain, just floor the difference!
          precalc_ndff += (fo_data.sigma - fi_data.sigma - 1)/n_phases + (fo_data.type == SA_GATE); //extra DFF before SA gate
        }
        return;
      }

      DEBUG_PRINT("\t[DFF] Non-straight chain: {}[{}] -> {}[{}]\n", GATE_TYPE.at(fi_data.type), (int)fi_data.sigma, GATE_TYPE.at(fo_data.type), (int)fo_data.sigma);
      
      // DFF locations from the earliest to the latest
      std::vector<uint64_t> dff_hashes;
      DEBUG_PRINT( "\tAdding new DFFs [reg size = {}] between phase {} and phase {}\n", DFF_REG.variables.size(), earliest_sigma, latest_sigma );

      for (glob_phase_t sigma = earliest_sigma; sigma <= latest_sigma; ++sigma)
      {
        uint64_t new_hash = DFF_REG.add(fi_node, fo_node, sigma);
        dff_hashes.push_back(new_hash);
        DEBUG_PRINT("\tAdded new DFFs at phase {} [reg size = {}]\n", sigma, DFF_REG.variables.size());
      }
      DEBUG_PRINT("\tConnecting new DFFs\n");
      for (auto i = 1u; i < dff_hashes.size(); ++i)
      {
        DFF_var & dff = DFF_REG.at( dff_hashes[i] );
        dff.parent_hashes.emplace(dff_hashes[i-1]);
      }

      if ( is_t1 )
      {
        DEBUG_PRINT("\tNode {} is a T1 gate, ", fo_node );
        DFF_closest_to_t1[fanin_ind] = dff_hashes.back();
        DEBUG_PRINT("\t{} is the DFF variable closest to this T1, on the path of the {}th input.\n", DFF_REG.str( dff_hashes.back() ) , ++fanin_ind );
      }

      // ensure that the SA gate is placed at least one stage after the fanin 
      if (fo_data.type == SA_GATE)
      {
        //(there has to be at least one DFF)
        assert( !dff_hashes.empty() );
        // The last DFF in the chain is required
        required_SA_DFFs.push_back(dff_hashes.back());
      }

      // if there are DFFs, the earliest_hash is the first hash in the chain
      // otherwise, it is the earliest_hash of the previous chain
      uint64_t earliest_hash = (dff_hashes.empty()) ? earliest_child_hash : dff_hashes.front();
      // if the node is internal, connect with the fanout phase
      if (fo_data.type == AA_GATE && !dff_hashes.empty() && earliest_hash != 0 && earliest_child_hash != 0)
      {
        DFF_var & child_dff = DFF_REG.at( earliest_child_hash );
        DEBUG_PRINT("\tPrior node is {}[{}]\n", child_dff.str(), (int)child_dff.sigma); 
        // assert(child_dff.fanin == fo_id);
        child_dff.parent_hashes.emplace( dff_hashes.back() );
      }
      if (fi_data.type == AA_GATE)
      {
        stack.emplace_back( fi_node, earliest_hash );
        DEBUG_PRINT("\tEmplacing {}({})[{}], {}\n", GATE_TYPE.at(fi_data.type), fi_node, (int)fi_data.sigma, (earliest_hash!=0)?DFF_REG.at(earliest_hash).str():"");
      }
    });

    if ( is_t1 )
    {
      DFF_closest_to_t1s.emplace( fo_node, DFF_closest_to_t1 );
    }
  }
  return std::make_tuple(DFF_REG, precalc_ndff, required_SA_DFFs);
}
