#include "L1Trigger/L1TMuonEndCap/interface/DebugTools.h"

namespace emtf {

  void dump_fw_raw_input(const l1t::EMTFHitCollection& out_hits, const l1t::EMTFTrackCollection& out_tracks) {
    // from interface/Common.h
    typedef L1TMuon::TriggerPrimitive TriggerPrimitive;
    constexpr int MIN_ENDCAP = 1;
    constexpr int MAX_ENDCAP = 2;
    constexpr int MIN_TRIGSECTOR = 1;
    constexpr int MAX_TRIGSECTOR = 6;

    // Input variables
    //
    // Field  | Value for  | Value for       | Value for
    // number | CSC        | RPC             | GE1/1
    // -------|------------|-----------------|----------------
    // 0      | bx_jitter  | 0               | 0
    // 1      | endcap     | endcap          | endcap
    // 2      | sector     | sector          | sector
    // 3      | subsector  | 0               | 0
    // 4      | station    | 10 deg chamber  | 10 deg chamber
    //        |            | (0..6)          | (0..6)
    // 5      | valid (=1) | valid (=2)      | valid (=3)
    // 6      | quality    | 0               | 0
    // 7      | pattern    | 0               | cluster size
    // 8      | wiregroup  | converted theta | partition
    // 9      | cscid      | station/ring    | layer
    //        |            | (0..5)          | (0..1)
    // 10     | bend       | 0               | 0
    // 11     | halfstrip  | converted phi   | pad

    // Delay w.r.t CSC
    //
    // Delay | RPC   | GE1/1
    // ------|-------|-------
    //       | +6    | -3

    for (int endcap = MIN_ENDCAP; endcap <= MAX_ENDCAP; ++endcap) {
      for (int sector = MIN_TRIGSECTOR; sector <= MAX_TRIGSECTOR; ++sector) {
        const int es = (endcap - MIN_ENDCAP) * (MAX_TRIGSECTOR - MIN_TRIGSECTOR + 1) + (sector - MIN_TRIGSECTOR);

        // _____________________________________________________________________
        // This prints the hits as raw text input to the firmware simulator
        // "12345" is the BX separator

        std::cout << "==== Endcap " << endcap << " Sector " << sector << " Hits ====" << std::endl;
        std::cout << "bx e s ss st vf ql cp wg id bd hs" << std::endl;

        bool valid_sector = std::any_of(out_hits.begin(), out_hits.end(), [&es](const auto& h) {
          return (h.Sector_idx() == es && h.Subsystem() == TriggerPrimitive::kCSC);
        });

        for (int ibx = -3 - 8; (ibx <= 3 + 11) && valid_sector; ++ibx) {
          for (const auto& h : out_hits) {
            if (h.Subsystem() == TriggerPrimitive::kCSC) {
              if (h.Sector_idx() != es)
                continue;
              if (!(-3 <= h.BX() && h.BX() <= 3 && h.BX() == ibx))
                continue;

              std::array<int, 12> values;
              values.fill(0);

              values[0] = 1;                                                                   // bx_jitter
              values[1] = (h.Endcap() == 1) ? 1 : 2;                                           // endcap
              values[2] = h.PC_sector();                                                       // sector
              values[3] = h.Subsector();                                                       // subsector
              values[4] = (h.PC_station() == 0 && h.Subsector() == 1) ? 1 : h.PC_station();    // station
              values[5] = 1;                                                                   // valid
              values[6] = h.Quality();                                                         // quality
              values[7] = h.Pattern();                                                         // pattern
              values[8] = h.Wire();                                                            // wiregroup
              values[9] = h.PC_chamber() + 1;                                                  // cscid
              values[10] = h.Bend();                                                           // bend
              values[11] = (h.Station() == 1 && h.Ring() == 4) ? h.Strip() + 128 : h.Strip();  // halfstrip

              for (unsigned i = 0; i < values.size(); i++) {
                if (i < values.size() - 1) {
                  std::cout << values[i] << " ";
                } else {
                  std::cout << values[i] << std::endl;
                }
              }

            } else if (h.Subsystem() == TriggerPrimitive::kRPC) {
              if (h.Sector_idx() != es)
                continue;
              if ((h.Station() == 3 || h.Station() == 4) && (h.Ring() == 1))
                continue;  // Skip iRPC
              if (!(-3 <= h.BX() && h.BX() <= 3 && (h.BX() + 6) == ibx))
                continue;  // RPC hits should be supplied 6 BX later relative to CSC hits

              // Assign RPC link index. Code taken from src/PrimitiveSelection.cc
              int rpc_sub = -1;
              int rpc_chm = -1;
              if (!h.Neighbor()) {
                rpc_sub = ((h.Subsector_RPC() + 3) % 6);
              } else {
                rpc_sub = 6;
              }
              if (h.Station() <= 2) {
                rpc_chm = (h.Station() - 1);
              } else {
                rpc_chm = 2 + (h.Station() - 3) * 2 + (h.Ring() - 2);
              }
              assert(rpc_sub != -1 && rpc_chm != -1);

              std::array<int, 12> values;
              values.fill(0);

              values[0] = 0;
              values[1] = (h.Endcap() == 1) ? 1 : 2;  // endcap
              values[2] = h.PC_sector();              // sector
              values[3] = 0;
              values[4] = rpc_sub;  // 10 deg chamber (0..6)
              values[5] = 2;        // valid (=2)
              values[6] = 0;
              values[7] = 0;
              values[8] = (h.Theta_fp() >> 2);  // converted theta
              values[9] = rpc_chm;              // station/ring (0..5)
              values[10] = 0;
              values[11] = (h.Phi_fp() >> 2);  // converted phi

              for (unsigned i = 0; i < values.size(); i++) {
                if (i < values.size() - 1) {
                  std::cout << values[i] << " ";
                } else {
                  std::cout << values[i] << std::endl;
                }
              }

            } else if (h.Subsystem() == TriggerPrimitive::kGEM) {
              if (h.Sector_idx() != es)
                continue;
              if (h.Station() == 2)
                continue;  // Skip GE2/1
              if (!(-3 <= h.BX() && h.BX() <= 3 && (h.BX() - 3) == ibx))
                continue;  // GEM hits should be supplied 3 BX earlier relative to CSC hits

              int gem_sub = -1;
              int gem_lay = -1;
              if (!h.Neighbor()) {
                gem_sub = h.PC_chamber();
              } else {
                gem_sub = 6;
              }
              const GEMDetId gemDetId = h.GEM_DetId();
              gem_lay = gemDetId.layer() - 1;
              assert(gem_sub != -1 && gem_lay != -1);

              std::array<int, 12> values;
              values.fill(0);

              values[0] = 0;
              values[1] = (h.Endcap() == 1) ? 1 : 2;  // endcap
              values[2] = h.PC_sector();              // sector
              values[3] = 0;
              values[4] = gem_sub;  // 10 deg chamber (0..6)
              values[5] = 3;        // valid (=3)
              values[6] = 0;
              values[7] = h.Quality();  // cluster size
              values[8] = h.Roll();     // partition
              values[9] = gem_lay;      // layer (0..1)
              values[10] = 0;
              values[11] = h.Strip();  // pad

              for (unsigned i = 0; i < values.size(); i++) {
                if (i < values.size() - 1) {
                  std::cout << values[i] << " ";
                } else {
                  std::cout << values[i] << std::endl;
                }
              }

            }  // end else-if
          }    // end loop over hits

          std::cout << "12345" << std::endl;  // BX separator

        }  // end loop over bx

        // _____________________________________________________________________
        // This prints the tracks as raw text output from the firmware simulator

        std::cout << "==== Endcap " << endcap << " Sector " << sector << " Tracks ====" << std::endl;
        std::cout << "bx e s a mo et ph cr q pt" << std::endl;

        for (const auto& t : out_tracks) {
          if (t.Sector_idx() != es)
            continue;

          std::cout << t.BX() << " " << (t.Endcap() == 1 ? 1 : 2) << " " << t.Sector() << " " << t.PtLUT().address
                    << " " << t.Mode() << " " << (t.GMT_eta() >= 0 ? t.GMT_eta() : t.GMT_eta() + 512) << " "
                    << t.GMT_phi() << " " << t.GMT_charge() << " " << t.GMT_quality() << " " << t.Pt() << std::endl;
        }  // end loop over tracks

      }  // end loop over sector
    }    // end loop over endcap
  }      // end function

}  // namespace emtf
