#include <cmath>
#include <memory>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
//#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"
//#include "L1Trigger/CSCCommonTrigger/interface/CSCConstants.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"

#include "helper.h"  // deltaPhiInDegrees

class MakeCoordLUTGE11 : public edm::EDAnalyzer {
public:
  explicit MakeCoordLUTGE11(const edm::ParameterSet&);
  virtual ~MakeCoordLUTGE11();

private:
  //virtual void beginJob();
  //virtual void endJob();

  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void endRun(const edm::Run&, const edm::EventSetup&);

  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  // Generate LUTs
  void generateLUTs();
  void generateLUTs_init();
  void generateLUTs_run();
  void generateLUTs_final();

  // Validate LUTs
  void validateLUTs();

  // Write LUT files
  void writeFiles();

  // Construct GEMDetId
  GEMDetId getGEMDetId(int endcap, int sector, int subsector, int station, int cscid, int layer, int roll) const;

  // Get global phi in degrees
  double getGlobalPhi(int endcap, int sector, int subsector, int station, int cscid, int layer, int roll, int pad) const;

  // Get global theta in degrees
  double getGlobalTheta(
      int endcap, int sector, int subsector, int station, int cscid, int layer, int roll, int pad) const;

  // Get sector phi in degrees
  double getSectorPhi(
      int endcap, int sector, int subsector, int station, int cscid, int layer, int roll, int pad, bool is_neighbor)
      const;

private:
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> me0_geom_token_;

  const edm::ParameterSet config_;

  int verbose_;

  std::string outdir_;

  bool please_validate_;

  int verbose_sector_;

  // Event setup
  const GEMGeometry* theGEMGeometry_;

  // Constants
  // [sector_12][station_1][chamber_7][layer_2]
  int ph_init[12][1][7][2];
  int ph_cover[12][1][7][2];
  int ph_reverse[12][1][7][2];

  // [sector_12][station_1][chamber_7][layer_2][roll_8]
  int th_mem[12][1][7][2][8];
};

// _____________________________________________________________________________
#define MIN_ENDCAP 1
#define MAX_ENDCAP 2
#define MIN_TRIGSECTOR 1
#define MAX_TRIGSECTOR 6

#define LOWER_THETA 8.5
#define UPPER_THETA 45.0

#define RAD_TO_DEG (180. / M_PI)

MakeCoordLUTGE11::MakeCoordLUTGE11(const edm::ParameterSet& iConfig)
    : me0_geom_token_(esConsumes<GEMGeometry, MuonGeometryRecord, edm::Transition::BeginRun>()),
      config_(iConfig),
      verbose_(iConfig.getUntrackedParameter<int>("verbosity")),
      outdir_(iConfig.getParameter<std::string>("outdir")),
      please_validate_(iConfig.getParameter<bool>("please_validate")),
      verbose_sector_(2) {
  // Zero multi-dimensional arrays
  memset(ph_init, 0, sizeof(ph_init));
  memset(ph_cover, 0, sizeof(ph_cover));
  memset(ph_reverse, 0, sizeof(ph_reverse));
  memset(th_mem, 0, sizeof(th_mem));
}

MakeCoordLUTGE11::~MakeCoordLUTGE11() {}

void MakeCoordLUTGE11::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
  // Geometry setup
  edm::ESHandle<GEMGeometry> gemGeometryHandle = iSetup.getHandle(me0_geom_token_);

  if (!gemGeometryHandle.isValid()) {
    std::cout << "ERROR: Unable to get MuonGeometryRecord!" << std::endl;
  } else {
    theGEMGeometry_ = gemGeometryHandle.product();
  }
}

void MakeCoordLUTGE11::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {}

void MakeCoordLUTGE11::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  generateLUTs();
  if (please_validate_)
    validateLUTs();
  writeFiles();
}

// _____________________________________________________________________________
void MakeCoordLUTGE11::generateLUTs() {
  generateLUTs_init();
  generateLUTs_run();
  generateLUTs_final();
}

void MakeCoordLUTGE11::generateLUTs_init() {
  for (const auto& etaPart : theGEMGeometry_->etaPartitions()) {
    //if (!etaPart->isGE11())
    //  continue;

    // Calculate phi-pitch from pitch, using parameters from GEMEtaPartitionSpecs
    auto pars = static_cast<const GEMEtaPartitionSpecs*>(&etaPart->type())->parameters();
    float b = pars.at(0);
    float B = pars.at(1);
    float h = pars.at(2);
    float r0 = h * (B + b) / (B - b);

    // Find z coordinate
    const LocalPoint& lp = etaPart->centreOfPad(etaPart->npads() / 2);
    const GlobalPoint& gp = etaPart->surface().toGlobal(lp);

    std::cout << "::generateLUTs_init()" << etaPart->id() << " " << etaPart->npads() << " " << etaPart->padPitch()
              << " " << std::setprecision(9) << (etaPart->padPitch() * RAD_TO_DEG) / r0 << " " << gp.z() << std::endl;
  }
  return;
}

// ph_init_hard values hardcoded in verilog zones.v
// These are in reduced precision, with negative offset relative to true strip 0.
// Negative offset must be large enough to accommodate chamber movement.
// [station_5][chamber_16]
// ME1 chambers 13,14,15,16 are neighbor sector chambers 3,6,9,12
// ME2 chambers 10,11 are neighbor sector chambers 3,9
static const int ph_init_hard[5][16] = {{39, 57, 76, 39, 58, 76, 41, 60, 79, 39, 57, 76, 21, 21, 23, 21},
                                        {95, 114, 132, 95, 114, 133, 98, 116, 135, 95, 114, 132, 0, 0, 0, 0},
                                        {38, 76, 113, 39, 58, 76, 95, 114, 132, 1, 21, 0, 0, 0, 0, 0},
                                        {38, 76, 113, 39, 58, 76, 95, 114, 132, 1, 21, 0, 0, 0, 0, 0},
                                        {38, 76, 113, 38, 57, 76, 95, 113, 132, 1, 20, 0, 0, 0, 0, 0}};

void MakeCoordLUTGE11::generateLUTs_run() {
  // theta_scale = 0.28515625 (7 bits to encode 128 values)
  constexpr double theta_scale = (UPPER_THETA - LOWER_THETA) / 128;
  // nominal_pitch = 0.133333 (ME2/2 strip pitch. 10-degree chamber, 80 strips - 5 overlap strips)
  constexpr double nominal_pitch = 10. / 75.;
  // pad_pitch ~ 0.0520833 (GE1/1 trigger pad pitch. 10-degree chamber, 192 pads)
  //constexpr double pad_pitch = 0.0529872;

  for (int endcap = MIN_ENDCAP; endcap <= MAX_ENDCAP; ++endcap) {
    for (int sector = MIN_TRIGSECTOR; sector <= MAX_TRIGSECTOR; ++sector) {
      for (int station = 1; station <= 1; ++station) {
        for (int subsector = 0; subsector <= 0; ++subsector) {  // unused
          for (int chamber = 1; chamber <= 7; ++chamber) {
            for (int layer = 1; layer <= 2; ++layer) {
              const bool is_neighbor = (chamber == 7) ? true : false;

              // Set 'real' CSCID, sector, subsector
              int rcscid = (chamber <= 3) ? chamber : chamber - 3;
              int rsector = sector;
              int rsubsector = (chamber <= 3) ? 1 : 2;

              if (is_neighbor) {
                rcscid = 3;
                rsector = (sector == 1) ? 6 : sector - 1;
                rsubsector = 2;
              }

              // Set firmware endcap/sector, station, chamber
              const int es = (endcap - 1) * 6 + (sector - 1);
              const int st = (station - 1);
              const int ch = (chamber - 1);
              const int ly = (layer - 1);
              assert(es < 12 && st < 1 && ch < 7 && ly < 2);

              // Specify maxRoll, maxPad
              const int maxRoll = 8;
              const int maxPad = 192;

              // find phi
              double fph_first =
                  getSectorPhi(endcap, rsector, rsubsector, station, rcscid, layer, maxRoll / 2, 0, is_neighbor);
              double fph_last = getSectorPhi(
                  endcap, rsector, rsubsector, station, rcscid, layer, maxRoll / 2, maxPad - 1, is_neighbor);
              assert(fph_first > 0. && fph_last > 0.);

              // set ph_init (phi coordinate of pad 0 in chamber)
              int my_ph_init_full = static_cast<int>(std::round(fph_first / (nominal_pitch / 8.)));  // 1/8-strip pitch
              ph_init[es][st][ch][ly] = my_ph_init_full;

              // set ph_cover (phi coordinate of pad maxPad-1 in chamber)
              int my_ph_cover_full = static_cast<int>(std::round(fph_last / (nominal_pitch / 8.)));  // 1/8-strip pitch
              ph_cover[es][st][ch][ly] = my_ph_cover_full;

              // set ph_reverse (0 if phi grows when pad # grows, 1 otherwise)
              int my_ph_reverse = deltaPhiInDegrees(fph_last, fph_first) > 0. ? 0 : 1;
              ph_reverse[es][st][ch][ly] = my_ph_reverse;

              if (verbose_ > 0 && sector == verbose_sector_) {
                double fph_first_glb =
                    getGlobalPhi(endcap, rsector, rsubsector, station, rcscid, layer, maxRoll / 2, 0);
                double fph_last_glb =
                    getGlobalPhi(endcap, rsector, rsubsector, station, rcscid, layer, maxRoll / 2, maxPad - 1);
                std::cout << "::generateLUTs_run()"
                          << " -- end " << endcap << " sec " << sector << " st " << st << " ch " << ch << " ly " << ly
                          << " -- fph_first_glb: " << fph_first_glb << " fph_last_glb: " << fph_last_glb
                          << " fph_first: " << fph_first << " fph_last: " << fph_last << std::endl;
              }

              // Loop over roll numbers
              for (int roll = 0; roll < maxRoll; ++roll) {
                // find theta for each roll
                // roll number starts from 1 in GEMDetId
                double fth_first = getGlobalTheta(endcap, rsector, rsubsector, station, rcscid, layer, roll + 1, 0);
                double fth_last =
                    getGlobalTheta(endcap, rsector, rsubsector, station, rcscid, layer, roll + 1, maxPad - 1);
                double fth = 0.5 * (fth_first + fth_last);
                assert(fth_first > 0. && fth_last > 0.);

                // set th_mem
                int my_th_mem = static_cast<int>(std::round((fth - LOWER_THETA) / theta_scale));
                th_mem[es][st][ch][ly][roll] = my_th_mem;
              }  // end loop over roll
            }    // end loop over layer
          }      // end loop over chamber
        }        // end loop over subsector
      }          // end loop over station
    }            // end loop over sector
  }              // end loop over endcap
  return;
}

void MakeCoordLUTGE11::generateLUTs_final() {
  for (int es = 0; es < 12; ++es) {
    for (int st = 0; st < 1; ++st) {
      for (int ch = 0; ch < 7; ++ch) {
        for (int ly = 0; ly < 2; ++ly) {
          assert(es < 12 && st < 1 && ch < 7 && ly < 2);

          // ph_init_hard is used to calculate zone_hit in the firmware
          // The following conditions must be satisfied to ensure that the logic
          //     ph_hit = ((fph + (1 << 4)) >> 5) - ph_init_hard;
          //     zone_hit = ph_hit + ph_init_hard;
          // yield
          //     zone_hit = ((fph + (1 << 4)) >> 5);

          int my_ph_init_hard = 0;
          if (ch < 3) {
            my_ph_init_hard = ph_init_hard[0][ch];
          } else if (ch < 6) {
            my_ph_init_hard = ph_init_hard[1][ch - 3];
          } else if (ch < 7) {
            my_ph_init_hard = ph_init_hard[0][12];
          }

          assert(((ph_init[es][st][ch][ly] + (1 << 4)) >> 5) >= my_ph_init_hard);
          assert(((ph_cover[es][st][ch][ly] + (1 << 4)) >> 5) >= my_ph_init_hard);

          // The following conditions must be satisfied to ensure that fph fits within 13 bits
          assert((ph_init[es][st][ch][ly] + (1 << 4)) < (1 << 13));
          assert((ph_cover[es][st][ch][ly] + (1 << 4)) < (1 << 13));
        }  // end loop over ly
      }    // end loop over ch
    }      // end loop over st
  }        // end loop over es
  return;
}

// Compare simulated (with floating-point) vs emulated (fixed-point) phi and theta coordinates
void MakeCoordLUTGE11::validateLUTs() {
  std::stringstream filename;

  filename << outdir_ << "/"
           << "validate.root";
  TFile* tfile = TFile::Open(filename.str().c_str(), "RECREATE");
  filename.str("");
  filename.clear();

  // Create TTree
  int lut_id = 0;
  int es = 0;
  int st = 0;
  int ch = 0;
  int ly = 0;
  //
  int endcap = 0;
  int station = 0;
  int sector = 0;
  int subsector = 0;
  int ring = 0;
  int chamber = 0;
  int CSC_ID = 0;
  int layer = 0;
  //
  int strip = 0;  // it is half-strip, despite the name
  int wire = 0;   // it is wiregroup, despite the name
  int fph_int = 0;
  int fth_int = 0;
  double fph_emu = 0.;  // in degrees
  double fth_emu = 0.;  // in degrees
  double fph_sim = 0.;  // in degrees
  double fth_sim = 0.;  // in degrees

  TTree* ttree = new TTree("tree", "tree");
  ttree->Branch("lut_id", &lut_id);
  ttree->Branch("es", &es);
  ttree->Branch("st", &st);
  ttree->Branch("ch", &ch);
  ttree->Branch("ly", &ly);
  //
  ttree->Branch("endcap", &endcap);
  ttree->Branch("station", &station);
  ttree->Branch("sector", &sector);
  ttree->Branch("subsector", &subsector);
  ttree->Branch("ring", &ring);
  ttree->Branch("chamber", &chamber);
  ttree->Branch("CSC_ID", &CSC_ID);
  ttree->Branch("layer", &layer);
  //
  ttree->Branch("strip", &strip);
  ttree->Branch("wire", &wire);
  ttree->Branch("fph_int", &fph_int);
  ttree->Branch("fth_int", &fth_int);
  ttree->Branch("fph_emu", &fph_emu);
  ttree->Branch("fth_emu", &fth_emu);
  ttree->Branch("fph_sim", &fph_sim);
  ttree->Branch("fth_sim", &fth_sim);

  for (es = 0; es < 12; ++es) {
    endcap = (es / 6) + 1;
    sector = (es % 6) + 1;

    for (st = 0; st < 1; ++st) {
      for (ch = 0; ch < 7; ++ch) {
        for (ly = 0; ly < 2; ++ly) {
          assert(es < 12 && st < 1 && ch < 7 && ly < 2);

          subsector = (ch < 3) ? 1 : 2;
          station = 1;
          ring = 1;
          chamber = ch + 1;
          layer = ly + 1;

          const bool is_neighbor = (chamber == 7) ? true : false;

          // Set 'real' CSCID, sector, subsector
          int rcscid = (chamber <= 3) ? chamber : chamber - 3;
          int rsector = sector;
          int rsubsector = (chamber <= 3) ? 1 : 2;

          if (is_neighbor) {
            rcscid = 3;
            rsector = (sector == 1) ? 6 : sector - 1;
            rsubsector = 2;
          }

          CSC_ID = rcscid;

          // Specify maxRoll, maxPad
          const int maxRoll = 8;
          const int maxPad = 192;

          // _______________________________________________________________________
          // Adapt logic from PrimitiveConversion

          // GE1/1 trigger pad pitch
          // pad_pitch_multiplier = round(pad_pitch / (nominal_pitch / 8) * 1024)
          const int pad_pitch_multiplier = 3256;  // 'factor' in PrimitiveConversion

          for (int roll = 0; roll < maxRoll; ++roll) {
            assert(es < 12 && st < 1 && ch < 7 && ly < 2 && roll < 8);

            for (int pad = 0; pad < maxPad; ++pad) {
              strip = pad;
              wire = roll;

              // ___________________________________________________________________
              // phi conversion

              int fph = ph_init[es][st][ch][ly];
              int fph_step_sign = (ph_reverse[es][st][ch][ly] == 0) ? 1 : -1;  // 'ph_tmp_sign' in PrimitiveConversion
              int fph_step = (pad * pad_pitch_multiplier) >> 10;               // 'ph_tmp' in PrimitiveConversion

              fph = fph + fph_step_sign * fph_step;

              // ___________________________________________________________________
              // theta conversion

              int th = th_mem[es][st][ch][ly][roll];

              // Protect against invalid value
              th = (th == 0) ? 1 : th;

              // ___________________________________________________________________
              // Finally

              // emulated phi and theta coordinates from fixed-point operations
              fph_int = fph;
              fph_emu = static_cast<double>(fph_int);
              fph_emu = fph_emu / 60.;
              fph_emu = fph_emu - 22. + 15. + (60. * (sector - 1));
              fph_emu = deltaPhiInDegrees(fph_emu, 0.);  // reduce to [-180,180]

              fth_int = th;
              fth_emu = static_cast<double>(fth_int);
              fth_emu = (fth_emu * (45.0 - 8.5) / 128. + 8.5);

              // simulated phi and theta coordinates from floating-point operations
              fph_sim = getGlobalPhi(endcap, rsector, rsubsector, station, rcscid, layer, roll + 1, pad);
              fth_sim = getGlobalTheta(endcap, rsector, rsubsector, station, rcscid, layer, roll + 1, pad);

              ttree->Fill();

              if (verbose_ > 1 && sector == verbose_sector_) {
                std::cout << "::validateLUTs()"
                          << " -- end " << endcap << " sec " << sector << " st " << st << " ch " << ch << " ly " << ly
                          << " wire " << wire << " strip " << strip << " -- fph_int: " << fph_int
                          << " fph_emu: " << fph_emu << " fph_sim: " << fph_sim << " -- fth_int: " << fth_int
                          << " fth_emu: " << fth_emu << " fth_sim: " << fth_sim << std::endl;
              }
            }  // end loop over pad
          }    // end loop over roll
        }      // end loop over ly
      }        // end loop over ch
    }          // end loop over st
  }            // end loop over es

  ttree->Write();
  tfile->Close();
  return;
}

// produce the LUT text files
void MakeCoordLUTGE11::writeFiles() {
  int num_of_files = 0;

  std::stringstream filename;

  for (int es = 0; es < 12; ++es) {
    int endcap = (es / 6) + 1;
    int sector = (es % 6) + 1;

    // write files: ph_init, ph_reverse
    std::ofstream ph_init_fs;
    filename << outdir_ << "/"
             << "ph_init_endcap_" << endcap << "_sect_" << sector << ".lut";
    ph_init_fs.open(filename.str().c_str());
    filename.str("");
    filename.clear();

    std::ofstream ph_reverse_fs;
    filename << outdir_ << "/"
             << "ph_reverse_endcap_" << endcap << "_sect_" << sector << ".lut";
    ph_reverse_fs.open(filename.str().c_str());
    filename.str("");
    filename.clear();

    for (int st = 0; st < 1; ++st) {
      for (int ch = 0; ch < 7; ++ch) {
        for (int ly = 0; ly < 2; ++ly) {
          assert(es < 12 && st < 1 && ch < 7 && ly < 2);

          ph_init_fs << std::hex << ph_init[es][st][ch][ly] << std::endl;
          ph_reverse_fs << std::hex << ph_reverse[es][st][ch][ly] << std::endl;
        }  // end loop over ly
      }    // end loop over ch
    }      // end loop over st

    ph_init_fs.close();
    ++num_of_files;
    ph_reverse_fs.close();
    ++num_of_files;

    // write files: th_mem
    for (int st = 0; st < 1; ++st) {
      for (int ch = 0; ch < 7; ++ch) {
        for (int ly = 0; ly < 2; ++ly) {
          int chamber = ch + 1;
          int layer = ly + 1;

          std::ofstream th_mem_fs;
          filename << outdir_ << "/"
                   << "th_mem_endcap_" << endcap << "_sect_" << sector << "_ch_" << chamber << "_ly_" << layer
                   << ".lut";
          th_mem_fs.open(filename.str().c_str());
          filename.str("");
          filename.clear();

          const int maxRoll = 8;
          for (int roll = 0; roll < maxRoll; ++roll) {
            assert(es < 12 && st < 1 && ch < 7 && ly < 2 && roll < 8);

            th_mem_fs << std::hex << th_mem[es][st][ch][ly][roll] << std::endl;
          }  // end loop over roll

          th_mem_fs.close();
          ++num_of_files;
        }  // end loop over ly
      }    // end loop over ch
    }      // end loop over st
  }        // end loop over es

  std::cout << "[INFO] Generated " << num_of_files << " LUT files." << std::endl;

  // Expects 12 sectors x (14 th_mem + 2 ph_init/ph_reverse)
  assert(num_of_files == 12 * (14 + 2));
  return;
}

// _____________________________________________________________________________
GEMDetId MakeCoordLUTGE11::getGEMDetId(
    int endcap, int sector, int subsector, int station, int cscid, int layer, int roll) const {
  const int region = (endcap == 2) ? -1 : 1;
  const int chamber = CSCTriggerNumbering::chamberFromTriggerLabels(sector, subsector, station, cscid);
  const int ring = 1;
  GEMDetId gemDetId = GEMDetId(region, ring, station, layer, chamber, roll);
  return gemDetId;
}

double MakeCoordLUTGE11::getGlobalPhi(
    int endcap, int sector, int subsector, int station, int cscid, int layer, int roll, int pad) const {
  const GEMDetId gemDetId = getGEMDetId(endcap, sector, subsector, station, cscid, layer, roll);
  const GEMEtaPartition* etaPart = theGEMGeometry_->etaPartition(gemDetId);
  assert(etaPart != nullptr);  // failed to get GEM eta partition
  const LocalPoint& lp = etaPart->centreOfPad(pad);
  const GlobalPoint& gp = etaPart->surface().toGlobal(lp);

  double phi = gp.barePhi() * RAD_TO_DEG;
  return phi;
}

double MakeCoordLUTGE11::getGlobalTheta(
    int endcap, int sector, int subsector, int station, int cscid, int layer, int roll, int pad) const {
  const GEMDetId gemDetId = getGEMDetId(endcap, sector, subsector, station, cscid, layer, roll);
  const GEMEtaPartition* etaPart = theGEMGeometry_->etaPartition(gemDetId);
  assert(etaPart != nullptr);  // failed to get GEM eta partition
  const LocalPoint& lp = etaPart->centreOfPad(pad);
  const GlobalPoint& gp = etaPart->surface().toGlobal(lp);

  double theta = gp.theta() * RAD_TO_DEG;
  theta = (endcap == 2) ? (180. - theta) : theta;  // put theta in the range of 0 - 180 degrees for negative endcap
  return theta;
}

double MakeCoordLUTGE11::getSectorPhi(
    int endcap, int sector, int subsector, int station, int cscid, int layer, int roll, int pad, bool is_neighbor)
    const {
  double globalPhi = getGlobalPhi(endcap, sector, subsector, station, cscid, layer, roll, pad);

  // Sector starts at -22 deg from the sector boundary
  double sectorStartPhi = -22. + 15. + (60. * (sector - 1));
  if (is_neighbor) {
    // This chamber comes from the neighbor sector into the native sector
    // Use the native sector sectorStartPhi (+60 deg)
    sectorStartPhi += 60.;
  }
  if (sectorStartPhi > 180.)
    sectorStartPhi -= 360.;

  double res = deltaPhiInDegrees(globalPhi, sectorStartPhi);
  assert(res >= 0.);
  return res;
}

// DEFINE THIS AS A PLUG-IN
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MakeCoordLUTGE11);
