#include <iostream>
#include <iomanip>
#include <map>
#include "rinex_obs.h"
#include "cycle_slip.h"
using namespace std;

int main()
{
	cout << "=== Full Constellation CycleSlip Scan ===" << endl;
	auto epochs = parseRinexAllEpochs("../../data/BILL00USA_R_20240070000_01D_15S_MO.rnx");
	cout << "Parsed " << epochs.size() << " epochs" << endl;

	// Per-satellite tracking state
	map<int, double> lg_prev, mw_prev;
	map<int, bool> have_prev;
	map<int, int> lg_slip_count;
	map<int, int> mw_slip_count;
	map<int, int> lli_l1_count;
	map<int, int> lli_l2_count;
	map<int, int> epoch_count;
	map<int, int> prev_epoch;
	map<int, int> new_arc_count;   // NEW: tracks large gaps separately from slips

	// Short gap = still compare across it (brief dropout, could be a real slip).
	// Long gap  = treat as a fresh arc (satellite re-acquired); don't compare,
	// since the jump is expected and NOT indicative of a slip within one arc.
	const int MAX_SHORT_GAP = 4;   // epochs (4 * 15s = 60s at this file's rate)

	int epoch_idx = 0;
	for (const auto& ep : epochs)
	{
		epoch_idx++;
		for (const auto& rec : ep.gps_records)
		{
			if (!rec.has_phase_l1 ||
				!rec.has_phase_l2 ||
				!rec.has_c2)
				continue;

			int prn = rec.prn;
			epoch_count[prn]++;

			double lg = geometryFree(rec.phase_l1, rec.phase_l2);
			double mw = melbourneWubbena(
				rec.pseudorange,
				rec.pseudorange_c2,
				rec.phase_l1,
				rec.phase_l2
			);

			if (have_prev[prn])
			{
				int gap = epoch_idx - prev_epoch[prn];

				if (gap >= 1 && gap <= MAX_SHORT_GAP)
				{
					// Short interruption (or normal consecutive epoch) --
					// still compare; a jump here is a real candidate slip.
					double dLG = lg - lg_prev[prn];
					double dMW = mw - mw_prev[prn];

					if (std::abs(dLG) > 1.0)
						lg_slip_count[prn]++;
					if (std::abs(dMW) > 5.0)
						mw_slip_count[prn]++;
				}
				else if (gap > MAX_SHORT_GAP)
				{
					// Long gap -- satellite was out of view / lost for a
					// while. This is a NEW ARC, not a slip within one --
					// don't compare LG/MW across it, just note it happened.
					new_arc_count[prn]++;
				}
			}

			if (rec.lli_l1 & 1)
				lli_l1_count[prn]++;
			if (rec.lli_l2 & 1)
				lli_l2_count[prn]++;

			lg_prev[prn] = lg;
			mw_prev[prn] = mw;
			have_prev[prn] = true;
			prev_epoch[prn] = epoch_idx;
		}
	}

	cout << endl << "PRN | Epochs | NewArcs | LG slips | MW slips | LLI_L1 | LLI_L2" << endl;
	cout << "===+========+=========+==========+==========+========+=======" << endl;
	for (auto& [prn, n] : epoch_count) {
		cout << "G" << (prn < 10 ? "0" : "") << prn
			 << " | " << n
			 << " | " << new_arc_count[prn]
			 << " | " << lg_slip_count[prn]
			 << " | " << mw_slip_count[prn]
			 << " | " << lli_l1_count[prn]
			 << " | " << lli_l2_count[prn] << endl;
	}

	return 0;
}