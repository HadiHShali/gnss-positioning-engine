#pragma once
// rinex_obs.h : RINEX observation reader, version-aware (2.xx and 3.xx)
// Week16 Day2: Replaces the week 15 header.
// ObsRecord / OpsEpoch fields are unchanged, so every week15 app still compiles against it. 

#include <string>
#include <vector>

// One GPS satellite at one Epoch
struct ObsRecord {
	int prn = 0;
	double pseudorange = 0.0;         // L1 code (m):  rinex versopn3 -> C1W (fallback C1P/C1Y/C1C) | rinex version2 -> P1 (fallback C1)
	double pseudorange_c2 = 0.0;	  // L2 code (m):  rinex versopn3 -> C2W (fallback C2P/C2Y)     | rinex version2 -> P2
	bool has_c2 = false;

	double phase_l1 = 0.0;		     // L1 phase (cycles): rinex versopn3 L1C (fallback L1W/L1P) | rinex versopn2 L1
	double phase_l2 = 0.0;			// L2 phase (cycles): rinex versopn3 L2W (fallback L2P/L2Y) | rinex versopn2 L2
	bool has_phase_l1 = false;
	bool has_phase_l2 = false;
	int lli_l1 = 0;					// loss-of-lock indicator, bit 0 = possible slip
	int lli_l2 = 0;
};


// Within one epoch
struct ObsEpoch {
	int year = 0, month = 0, day = 0, hour = 0, minute = 0;
	double second = 0.0;
	double t_gps = 0.0; // seconds of GPS week
	int gps_week = 0;  // NEW: full GPS week number
	std::vector<ObsRecord> gps_records;
};

// header summary
struct RinexObsHeader {
	double version = 0.0;
	int major = 0;                 // 2, 3 or 4 (4.xx obs is read with the v3 layout)
	std::string marker_name;
	double approx_x = 0.0, approx_y = 0.0, approx_z = 0.0;   // APPROX POSITION XYZ (m)
	double interval = 0.0;         // INTERVAL (second); 0 if the header has none
	std::vector<std::string> gps_obs_types; // as written: "L1","P2" (v2) or "L1C","C2W" (v3)
	std::string code_p1, code_p2, code_l1, code_l2;   // observables actually used ("-" = none)

};

// Reads only line 1. Returns e.g. 2.11 or 3.04, or 0.0 if it is not a RINEX file.
// file_type (optional) receives 'O' (obs), 'N' (navigation), ...
double detectRinexVersion(const std::string& filename, char* file_type = nullptr);

// Parse every epoch; dispatches to the v2 or v3 reader from line1.
std::vector<ObsEpoch> parseRinexAllEpochs(const std::string& filename, RinexObsHeader* hdr = nullptr);