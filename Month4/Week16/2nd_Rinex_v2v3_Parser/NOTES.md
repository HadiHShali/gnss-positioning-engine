# Week 16 Day 2 — RINEX 2/3 version-aware parsers

**Date:** 2026-10-05 · **Data:** 2024-01-07 (DOY 007, GPS week 2296), 15 s
**Stations:** BILL (Week 15), LCHS (RTK base), NMKM (RTK rover)

## What changed

- `rinex_obs.cpp` rewritten: reads RINEX **2 and 3** obs (version from line 1).
  `ObsRecord`/`ObsEpoch` unchanged → old apps still compile. New: `RinexObsHeader`, `gps_week`.
- `satellites.cpp` patched: `parseRinexNav` and `parseKlobucharFromNav` read RINEX **2 and 3**;
  `findBestEphemeris` now skips unhealthy satellites (`sv_health != 0`) and records
  more than 2 h from Toe (4-hour fit interval).
- **Bug fixed:** old nav loop skipped 7 lines for every non-GPS record, but GLONASS/SBAS
  records are shorter → lost GPS records in mixed files (81 of 162 in a test; 0 from RINEX 2).
- `test_phase_parser_allSat` takes the obs file as an argument.
- New app `check_station`: one-station check of obs + nav before any differencing.
- RTK pair switched to **LCHS (base) / NMKM (rover)**.

## Files

| File | Origin |
|---|---|
| `include/rinex_obs.h`, `src/rinex_obs.cpp` | written today |
| `include/satellites.h`, `src/satellites.cpp` | from `Month3/Week14/3rd_SP3_SPP_Solver`, patched |
| `include/cycle_slip.h`, `apps/test_phase_parser_allSat.cpp` | from `Month4/Week15/3rd_Phase_Parser_AllSat` |
| `apps/check_station.cpp` | written today |

## Build and run (from this folder)

```bash
cmake -B build && cmake --build build --config Release
./build/Release/check_station.exe data/<obs> [data/<nav>]
./build/Release/test_phase_parser_allSat.exe data/<obs>
```

Repo is on OneDrive: "Permission denied" / "access denied" = file locked by sync → close
editors or pause sync. Quote full paths (they contain spaces).

## Data (EarthScope GAGE archive, login required)

| File | Version |
|---|---|
| `BILL00USA_R_..._15S_MO.rnx` / `bill0070.24o` | 3.04 / 2.11 obs |
| `bill0070.24n` | 2.11 nav (GPS only) |
| `LCHS00USA_R_..._15S_MO.rnx` + `_MN.rnx` | 3.04 obs + mixed nav |
| `NMKM00USA_R_..._15S_MO.rnx` + `_MN.rnx` | 3.04 obs + mixed nav |

RINEX 2: `/archive/gnss/rinex/obs/2024/007/bill0070.24d.Z` (take `.d.Z`; `.o.Z` expires).
```bash
gzip -d data/bill0070.24d.Z
/c/Dev/tools/RNXCMP_4.2.0_Windows_mingw_64bit/bin/CRX2RNX.exe data/bill0070.24d
```
Line 1 ending in `CRINEX VERS / TYPE` = still Hatanaka. First byte `0x1F` = still .gz/.Z.
The parsers print both messages.

## How the parsers work (key ideas)

- **Obs slot:** every value is 16 chars (14 value + LLI + SSI) in both versions.
  v3: one line per satellite → `line.substr(3)`. v2: 5 slots per line → join the
  satellite's lines, each **padded to 80**. Then slot `i` is at `16*i` for both.
- **Obs types:** v2 = one list for all systems, 2-char codes at `10+6i`, 9/line.
  v3 = one list per system, 3-char codes at `7+4i`, 13/line.
- **v2 epochs:** no `>`; satellite list on the epoch line (12/line); data lines have no ID —
  order is everything. 2-digit year: 80–99 → 19xx.
- **Epoch flags:** 2–5 skip `n` lines, 6 read and discard. Always read every line.
- **Codes used:** v3 C1W/C2W/L1C/L2W (fallbacks in P(Y) family; C1C last). v2 P1(C1)/P2/L1/L2.
  L2C codes never used as fallback.
- **Nav:** same 29 values in both; differences are the first-line layout
  (clock fields at 22 vs 23) and the orbit-line margin (**3 vs 4**).
  Exponent D/E is the producer's choice (teqc D, sbf2rin E).
- **GPS time:** days since 1980-01-06 → week = days/7, sow = (days%7)·86400 + time.

## check_station — why

Parse errors don't crash; they give plausible wrong numbers. The nav checks test physics:
A ≈ 26,560 km, e < 0.03, i ≈ 55°, **Ω̇ < 0**, nav week = obs week.
`stod()` reads the longest valid number and ignores the rest: one column too far right
keeps the digits but **drops the minus sign** — only the Ω̇ sign check catches that.

## Results (all checks OK)

| | BILL v2 | LCHS v3 | NMKM v3 |
|---|---|---|---|
| Codes | P1 P2 L1 L2 | C1W C2W L1C L2W | C1W C2W L1C L2W |
| Epochs / PRNs | 5760 / 30 | 5760 / 30 | 5760 / 30 |
| Ephemerides | 208 | 185 | 195 |
| Unhealthy | G01 | G01 | G01 |
| Records without usable eph | 0 / 59,653 | 0 / 51,067 | 0 / 52,896 |

**Baseline LCHS → NMKM: 10.79 km, azimuth 37°, Δh −15 m** (header XYZ, approximate).
Short: atmosphere and orbit errors mostly cancel. Use LCHS coordinates from PhD solutions.

## Findings

- **RINEX 2 "L2" is mixed:** BILL v2 L2 = v3 **L2L** (L2C) for G06, P(Y) for satellites
  without L2C. RINEX 2 can't tell. With the same app, v2 shows more LG slips than v3
  (G05 5→9, G07 0→5, G10 1→5, G30 1→6, G31 0→3) — likely L2C↔P(Y) switching
  (~1.96 m jump > 1 m threshold), unverified.
- **LLI flags depend on the converter** (same L1 data, different LLI_L1 counts).
- **Site effects:** BILL G20 and LCHS G09/G26/G29 are rough but normal at NMKM.
  Don't use G09/G26/G29 as reference satellite.
- **G01 unhealthy** (health 63) in all nav files and not in any obs file
  (`grep -c G01 data/bill0070.24o` = 0): the receiver/converter left it out.
- `bill0070.24o` has no INTERVAL line — measure it from the data.
- Single-station nav files have ~6–7 ephemerides per PRN (BRDC ~14).

## Corrections

- L2C vs L2W is **not** a quarter cycle in RINEX 3 (BILL: 8.011 cycles). Don't mix them
  because they're separate tracking loops.
- The Week 15 table IS correct (G20 `1|12|12`). The `2|11|11` numbers came from the Week 15
  folder's app, where `MAX_SHORT_GAP` had been changed to **3** with a stale comment
  ("4 * 15s = 60s"). Day 2 uses 4 and reproduces the Week 15 table exactly.
  Lesson: change a value's comment in the same edit.

## TODO results (done 2026-10-05)

- [x] Full diff old vs new parser on BILL v3 → identical (only the 4 threshold-dependent rows differed)
- [x] `MAX_SHORT_GAP`: Week 15 copy = 3, Day 2 = 4 → explains all table differences and the
      NewArcs "puzzle" (not a data effect). Reset the Week 15 copy to 4.
- [x] G01 in `bill0070.24o`? No (count 0).
- [x] `findBestEphemeris` age limit: 2 h (`MAX_AGE = 7200`, as RTKLIB). `check_station` now tests
      every record at its own time (a single mid-file time gave false alarms).
- [x] Klobuchar reader: RINEX 2 `ION ALPHA/BETA` + RINEX 3 `GPSA/GPSB`; arrays default to 0.
      Tested on RTKLIB `07590920.05n` (v2) and the LCHS/NMKM-style v3 header.

**Known limitations (no data to test):** RINEX 4 obs (read with v3 layout); RINEX 4 nav
(rejected); mid-file header changes (event flag 4) are skipped; GPS only.

## Next: RTK

Match epochs LCHS↔NMKM → code DD → phase DD (float) → LAMBDA in Week 17.
Reference satellite: high elevation, clean at both stations.
