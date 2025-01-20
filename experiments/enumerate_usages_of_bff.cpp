#include <cstdint>
#include <fmt/format.h>

#include <kitty/static_truth_table.hpp>
#include <kitty/npn.hpp>

typedef kitty::static_truth_table<3u> TT3;

struct n_config
{
	uint8_t truth_table;
  uint8_t phase;

  n_config() {}
  n_config( uint8_t truth_table_, uint32_t phase_ ) : truth_table( truth_table_ ), phase( phase_ ) {}
  n_config( n_config const& other ) : truth_table( other.truth_table ), phase( other.phase ) {}
};

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
  /* TTs of the 3 outputs */
  uint8_t sum_truth_table{};
  uint8_t carry_truth_table{};
  uint8_t cbar_truth_table{};

  // TODO : calculate cost, add mask for specific outputs
  T1_OUTPUTS() : in_phase(0u), has_carry_inverted(false), has_cbar_inverted(false), has_carry(false), has_cbar(false) {}

  // Explicit Constructor
  T1_OUTPUTS(const uint8_t in, const bool p_c, const bool p_cb, const bool h_c, const bool h_cb)
    : in_phase(in), has_carry_inverted(p_c), has_cbar_inverted(p_cb), has_carry(h_c), has_cbar(h_cb) {}

  /* Explicit constructor with output functions provided */
  T1_OUTPUTS( const uint8_t in, const bool p_c, const bool p_cb, const uint8_t tt_sum, const uint8_t tt_carry, const uint8_t tt_cbar )
  	: in_phase( in ), has_carry_inverted( p_c ), has_cbar_inverted( p_cb ), has_carry( true ), has_cbar( true ), sum_truth_table( tt_sum ), carry_truth_table( tt_carry ), cbar_truth_table( tt_cbar ) {}
    
  // // returns the TT of the sum output of the T1 cell
  // uint8_t sum_tt() const
  // {
  //   return (__builtin_popcount(in_phase) & 1) ? 0x69 : 0x96;
  // }
  // // returns the TT of the carry output of the T1 cell
  // uint8_t carry_tt(const bool output_phase = false) const
  // {
  //   TT3 tt_in;
  //   tt_in._bits = MAJ3;
  //   const uint32_t phase = in_phase | (output_phase << CUT_SIZE);
    
  //   TT3 tt_out = kitty::create_from_npn_config( std::make_tuple( tt_in, phase, perm ) );
  //   return static_cast<uint8_t>(tt_out._bits);
  // }
  // // returns the TT of the cbar output of the T1 cell
  // uint8_t cbar_tt(const bool output_phase = false) const
  // {
  //   TT3 tt_in;
  //   tt_in._bits = OR3;
  //   const uint32_t phase = in_phase | (output_phase << CUT_SIZE);
    
  //   TT3 tt_out = kitty::create_from_npn_config( std::make_tuple( tt_in, phase, perm ) );
  //   return static_cast<uint8_t>(tt_out._bits);
  // }
  // uint8_t cost() const
  // {
  //   return T1_COST // base cost (2x merger + T1 cell)
  //     + __builtin_popcount(in_phase) * (COST_MAP[fNOT] - COST_MAP[fDFF]) // cost of negating the inputs
  //     + COST_MAP[fDFF] * (has_carry + has_cbar)                          // cost of MAJ3 and OR3
  //     + COST_MAP[fNOT] * (has_carry_inverted + has_cbar_inverted)        // cost of ~(MAJ3) and ~(OR3)
  //     + COST_MAP[fSPL] * ((has_carry & has_carry_inverted) + (has_cbar & has_cbar_inverted)); //cost of splitters, if needed 
  // }
  void report() const
  {
    fmt::print("\t\tInput phases : {{ {}0,{}1,{}2 }}\n", ( ( in_phase & 1 ) ? "~" : " " ), ( ( in_phase >> 1 & 1 ) ? "~" : " " ), ( ( in_phase >> 2 & 1 ) ? "~" : " " ) );
    fmt::print("\t\tXOR3: {0:08b} ( 0x{0:02x} )\n",   sum_truth_table );
    fmt::print("\t\tMAJ3: {0:08b} ( 0x{0:02x} )\n", carry_truth_table );
    fmt::print("\t\t OR3: {0:08b} ( 0x{0:02x} )\n",  cbar_truth_table );
  }
};

int main()
{
	TT3 XOR3, MAJ3, OR3;
	XOR3._bits = 0x96;
	MAJ3._bits = 0xe8;
	OR3._bits  = 0xfe;

	std::vector<n_config> xor3_n_configs;
	std::vector<n_config> maj3_n_configs;
	std::vector<n_config>  or3_n_configs;
	auto add_xor3 = [&xor3_n_configs]( auto const& tt, auto const& phase ){ xor3_n_configs.emplace_back( tt._bits, phase ); };
	auto add_maj3 = [&maj3_n_configs]( auto const& tt, auto const& phase ){ maj3_n_configs.emplace_back( tt._bits, phase ); };
	auto add_or3  =  [&or3_n_configs]( auto const& tt, auto const& phase ){  or3_n_configs.emplace_back( tt._bits, phase ); };

	kitty::exact_n_enumeration( XOR3, add_xor3 );
  kitty::exact_n_enumeration( MAJ3, add_maj3 );
  kitty::exact_n_enumeration(  OR3,  add_or3 );

  std::vector<T1_OUTPUTS> BFF_usages;

  /* enumerate all the potential combinations */
  for ( auto const& xor3_config: xor3_n_configs )
  {
  	uint8_t phase_target{ xor3_config.phase };
  	for ( auto const& maj3_config: maj3_n_configs )
  	{
  		if ( maj3_config.phase != phase_target )
  		{
  			continue;
  		}
  		for ( auto const& or3_config: or3_n_configs )
  		{
  			if ( or3_config.phase != phase_target )
  			{
  				continue;
  			}

				BFF_usages.push_back( T1_OUTPUTS( phase_target, false, false, 
				                                  xor3_config.truth_table, 
				                                  maj3_config.truth_table, 
				                                  or3_config.truth_table ) );
  		}
  	}
  }

  /* print out all the potential usage of BFF */
  fmt::print( "There are {} cases in total. \n", BFF_usages.size() );
  uint32_t case_ctr{};
  for ( auto const& BFF_usages_each: BFF_usages )
  {
  	fmt::print( "\t[{}]: \n", ++case_ctr );
  	BFF_usages_each.report();
  	fmt::print( "\n" );
  }
}
\
	