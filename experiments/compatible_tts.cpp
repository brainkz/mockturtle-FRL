#include <iostream>
#include <bit>
#include <tuple>
#include <vector>
#include <unordered_map>

#include <mockturtle/algorithms/nodes.hpp>
#include <kitty/kitty.hpp>


constexpr uint8_t CUT_SIZE = 3u;
// constexpr uint64_t NUM_TTS = (1 << CUT_SIZE);
constexpr uint64_t NUM_OUT = (1 << CUT_SIZE);
constexpr auto COST_MAP = COSTS_SUNMAGNETICS;
constexpr uint8_t T1_COST = 9u + 2*COST_MAP[fCB];

constexpr uint8_t XOR3 = 0x96;
constexpr uint8_t MAJ3 = 0xe8;
constexpr uint8_t  OR3 = 0xfe;
const std::vector<uint8_t> perm { 0, 1, 2 };

struct T1_OUTPUTS
{
  // 3 bits to store input negations
  unsigned in_phase          : 3;
  // whether carry (maj3) output is used with inverter, feel free to rename if other naming is more convenient for you
  bool     has_carry_inverted: 1;
  // whether cbar (or3) output is used with inverter, feel free to rename if other naming is more convenient for you
  bool     has_cbar_inverted : 1;
  // whether carry (maj3) output is used with DFF
  bool     has_carry         : 1;
  // whether cbar (or3) output is used with DFF
  bool     has_cbar          : 1;

  // TODO : calculate cost, add mask for specific outputs
  T1_OUTPUTS() : in_phase(0u), has_carry_inverted(false), has_cbar_inverted(false), has_carry(false), has_cbar(false) {}

  // Explicit Constructor
  T1_OUTPUTS(const uint8_t in, const bool p_c, const bool p_cb, const bool h_c, const bool h_cb)
    : in_phase(in), has_carry_inverted(p_c), has_cbar_inverted(p_cb), has_carry(h_c), has_cbar(h_cb) {}
    
  // returns the TT of the sum output of the T1 cell
  uint8_t sum_tt() const
  {
    return (__builtin_popcount(in_phase) & 1) ? 0x69 : 0x96;
  }
  // returns the TT of the carry output of the T1 cell
  uint8_t carry_tt(const bool output_phase = false) const
  {
    TT3 tt_in;
    tt_in._bits = MAJ3;
    const uint32_t phase = in_phase | (output_phase << CUT_SIZE);
    
    TT3 tt_out = kitty::create_from_npn_config( std::make_tuple( tt_in, phase, perm ) );
    return static_cast<uint8_t>(tt_out._bits);
  }
  // returns the TT of the cbar output of the T1 cell
  uint8_t cbar_tt(const bool output_phase = false) const
  {
    TT3 tt_in;
    tt_in._bits = OR3;
    const uint32_t phase = in_phase | (output_phase << CUT_SIZE);
    
    TT3 tt_out = kitty::create_from_npn_config( std::make_tuple( tt_in, phase, perm ) );
    return static_cast<uint8_t>(tt_out._bits);
  }
  uint8_t cost() const
  {
    return T1_COST // base cost (2x merger + T1 cell)
      + __builtin_popcount(in_phase) * (COST_MAP[fNOT] - COST_MAP[fDFF]) // cost of negating the inputs
      + COST_MAP[fDFF] * (has_carry + has_cbar)                          // cost of MAJ3 and OR3
      + COST_MAP[fNOT] * (has_carry_inverted + has_cbar_inverted)        // cost of ~(MAJ3) and ~(OR3)
      + COST_MAP[fSPL] * ((has_carry & has_carry_inverted) + (has_cbar & has_cbar_inverted)); //cost of splitters, if needed 
  }
};

int main()
{
  std::array<T1_OUTPUTS, NUM_OUT> compatible_tts;
}

#if false 
  // OLD implementation, likely useless
  struct T1_OUTPUTS
  {
    uint8_t in_phase;
    TT3 sum_tt;
    TT3 carry_tt;
    TT3 cbar_tt;
    bool phase_carry; //false means no inversion
    bool phase_cbar; //false means no inversion
    bool has_carry;
    bool has_cbar;
    uint32_t cost; //each negative input phase increases the cost by (COST_MAP[fNOT] - COST_MAP[fDFF])

    // TODO : calculate cost, add mask for specific outputs
    T1_OUTPUTS() : in_phase(0u), phase_carry(false), phase_cbar(false), has_carry(false), has_cbar(false)
    {
      sum_tt._bits = 0x96;
      carry_tt._bits = 0xe8;
      cbar_tt._bits = 0xfe;
      cost = T1_COST;
    }

    // Explicit Constructor
    T1_OUTPUTS(const uint8_t in, const TT3& s, const TT3& c, const TT3& cb, const bool p_c, const bool p_cb, const bool h_c, const bool h_cb, const uint32_t cst)
      : in_phase(in), sum_tt(s), carry_tt(c), cbar_tt(cb), phase_carry(p_c), phase_cbar(p_cb), has_carry(h_c), has_cbar(h_cb), cost(cst) {}
      
    // Explicit Constructor
    T1_OUTPUTS(const uint8_t in, const TT3& s, const TT3& c, const TT3& cb, const bool p_c, const bool p_cb, const bool h_c, const bool h_cb)
      : in_phase(in), sum_tt(s), carry_tt(c), cbar_tt(cb), phase_carry(p_c), phase_cbar(p_cb), has_carry(h_c), has_cbar(h_cb) 
    {
      cost = T1_COST + __builtin_popcount(in) * (COST_MAP[fNOT] - COST_MAP[fDFF]) + COST_MAP[fDFF] * (has_carry & ~phase_carry + has_cbar & ~phase_cbar) + COST_MAP[fNOT] * (has_carry &  phase_carry + has_cbar &  phase_cbar);
    }
  };
#endif
