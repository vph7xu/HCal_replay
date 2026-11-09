// hcal_waveform_energy_study.C
//
// Stage-1 HCal waveform reconstruction/optimization study.
//
// Uses the existing SBS primary-cluster membership (sbs.hcal.clus_blk.id),
// but reconstructs the charge of each cluster block directly from
// sbs.hcal.samps using:
//   1) residual baseline estimate from waveform tails,
//   2) constant-fraction timing,
//   3) a fixed integration window relative to the CFD crossing.
//
// This lets us study CFD fraction and integration-window choices without
// changing the current cluster membership/threshold logic yet.
//
// Required branches:
//   sbs.hcal.samps
//   sbs.hcal.samps_idx
//   sbs.hcal.nsamps
//   sbs.hcal.samps_elemID
//   sbs.hcal.clus_blk.id
//   sbs.hcal.e
//
// Typical usage:
//   root [0] .L hcal_waveform_energy_study.C+
//   root [1] draw_hcal_waveform("replay.root",100,132)
//   root [2] hcal_waveform_energy_study("replay.root")
//   root [3] hcal_scan_windows("replay.root")
//
// Or directly from the shell:
//   root -l -q 'hcal_waveform_energy_study.C+("replay.root")'
//
// Notes:
//   * sbs.hcal.samps is documented as CALIBRATED ADC samples.
//   * This code still removes a small residual baseline measured from the
//     stored waveform tails. It does NOT subtract sbs.hcal.ped again.
//   * hcal_scan_windows() minimizes the width of Qcustom / sbs.hcal.e.
//     That is a reconstruction-consistency metric, not by itself a true
//     detector energy-resolution measurement. For a true resolution scan,
//     replace the reference with an independent expected-energy quantity.

#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

namespace HCalWF {

struct PulseResult {
  bool ok = false;
  int elemID = -1;
  int nsamp = 0;
  int peakBin = -1;

  double baseline = 0.0;     // residual baseline in stored/calibrated samples
  double noiseRMS = 0.0;
  double peak = 0.0;
  double snr = 0.0;

  double t10 = std::numeric_limits<double>::quiet_NaN(); // sample units
  double t20 = std::numeric_limits<double>::quiet_NaN();
  double tCFD = std::numeric_limits<double>::quiet_NaN();
  double t80 = std::numeric_limits<double>::quiet_NaN();
  double t90 = std::numeric_limits<double>::quiet_NaN();

  double rise1090 = std::numeric_limits<double>::quiet_NaN(); // sample units
  double rise2080 = std::numeric_limits<double>::quiet_NaN();

  int intLo = -1;
  int intHi = -1;
  double qFixed = 0.0;
  double qFull = 0.0;
};

static double mean(const std::vector<double>& v)
{
  if (v.empty()) return 0.0;
  return std::accumulate(v.begin(), v.end(), 0.0) / double(v.size());
}

static double rms_about(const std::vector<double>& v, double m)
{
  if (v.size() < 2) return 0.0;
  double s2 = 0.0;
  for (double x : v) {
    const double d = x - m;
    s2 += d*d;
  }
  return std::sqrt(s2 / double(v.size()-1));
}

// Estimate the residual baseline using whichever waveform tail has the lower
// average. For positive pulses this tends to choose the pulse-free side.
static void estimate_baseline(const std::vector<double>& w,
                              int nPed,
                              double& baseline,
                              double& noiseRMS)
{
  baseline = 0.0;
  noiseRMS = 0.0;
  const int n = int(w.size());
  if (n <= 0) return;

  int np = std::min(nPed, std::max(1, n/4));
  if (2*np > n) np = std::max(1, n/2);

  std::vector<double> first, last;
  first.reserve(np);
  last.reserve(np);
  for (int i = 0; i < np; i++) {
    first.push_back(w[i]);
    last.push_back(w[n-1-i]);
  }

  const double mf = mean(first);
  const double ml = mean(last);

  // Positive HCal pulses: choose the lower-mean tail, following the same
  // basic protection against a pulse contaminating one end of the waveform.
  const std::vector<double>& chosen = (mf <= ml) ? first : last;
  baseline = mean(chosen);
  noiseRMS = rms_about(chosen, baseline);
}

// Find the rising-edge crossing closest to the pulse peak at a requested
// fraction of the peak. Linear interpolation gives sub-sample timing.
static double crossing_time(const std::vector<double>& v,
                            int peakBin,
                            double fraction)
{
  if (peakBin <= 0 || peakBin >= int(v.size()))
    return std::numeric_limits<double>::quiet_NaN();

  const double level = fraction * v[peakBin];
  for (int i = peakBin-1; i >= 0; --i) {
    if (v[i] < level && v[i+1] >= level) {
      const double dy = v[i+1] - v[i];
      if (std::fabs(dy) < 1e-12) return double(i);
      return double(i) + (level - v[i]) / dy;
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

static PulseResult analyze_waveform(const TTreeReaderArray<Double_t>& samps,
                                    int start,
                                    int n,
                                    int elemID,
                                    double cfdFrac,
                                    int nPre,
                                    int nPost,
                                    int nPed,
                                    double minSNR = 0.0)
{
  PulseResult r;
  r.elemID = elemID;
  r.nsamp = n;

  if (start < 0 || n < 3 || start+n > int(samps.GetSize())) return r;
  if (!(cfdFrac > 0.0 && cfdFrac < 1.0)) return r;

  std::vector<double> w(n), v(n);
  for (int i = 0; i < n; ++i) w[i] = samps[start+i];

  estimate_baseline(w, nPed, r.baseline, r.noiseRMS);
  for (int i = 0; i < n; ++i) v[i] = w[i] - r.baseline;

  r.peakBin = int(std::max_element(v.begin(), v.end()) - v.begin());
  r.peak = v[r.peakBin];
  if (!(r.peak > 0.0)) return r;

  if (r.noiseRMS > 0.0) r.snr = r.peak / r.noiseRMS;
  else                  r.snr = std::numeric_limits<double>::infinity();

  if (minSNR > 0.0 && r.snr < minSNR) return r;

  r.t10  = crossing_time(v, r.peakBin, 0.10);
  r.t20  = crossing_time(v, r.peakBin, 0.20);
  r.tCFD = crossing_time(v, r.peakBin, cfdFrac);
  r.t80  = crossing_time(v, r.peakBin, 0.80);
  r.t90  = crossing_time(v, r.peakBin, 0.90);

  if (!std::isfinite(r.tCFD)) return r;

  if (std::isfinite(r.t10) && std::isfinite(r.t90)) r.rise1090 = r.t90-r.t10;
  if (std::isfinite(r.t20) && std::isfinite(r.t80)) r.rise2080 = r.t80-r.t20;

  // Anchor the integer integration window to the CFD crossing bin.
  const int anchor = int(std::floor(r.tCFD));
  r.intLo = std::max(0, anchor - nPre);
  r.intHi = std::min(n-1, anchor + nPost);

  for (int i = r.intLo; i <= r.intHi; ++i) r.qFixed += v[i];
  for (double x : v) r.qFull += x;

  r.ok = true;
  return r;
}

static int find_waveform_index(const TTreeReaderArray<Double_t>& elemID,
                               int wantedID)
{
  for (int i = 0; i < int(elemID.GetSize()); ++i) {
    if (int(std::llround(elemID[i])) == wantedID) return i;
  }
  return -1;
}

static TFile* open_file_and_tree(const char* infile, TTree*& T)
{
  T = nullptr;
  if (!infile || std::string(infile).empty()) {
    std::cerr << "ERROR: no input ROOT file was specified.\n";
    return nullptr;
  }

  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: could not open input ROOT file: " << infile << "\n";
    if (fin) { fin->Close(); delete fin; }
    return nullptr;
  }

  T = dynamic_cast<TTree*>(fin->Get("T"));
  if (!T) {
    std::cerr << "ERROR: tree T was not found in " << infile << "\n";
    fin->Close();
    delete fin;
    return nullptr;
  }
  return fin;
}

static bool has_required_branches(TTree* T)
{
  const char* req[] = {
    "sbs.hcal.samps",
    "sbs.hcal.samps_idx",
    "sbs.hcal.nsamps",
    "sbs.hcal.samps_elemID",
    "sbs.hcal.clus_blk.id",
    "sbs.hcal.e"
  };

  bool ok = true;
  for (const char* b : req) {
    if (!T->GetBranch(b)) {
      std::cerr << "ERROR: missing required branch: " << b << "\n";
      ok = false;
    }
  }
  return ok;
}

} // namespace HCalWF


// ---------------------------------------------------------------------------
// Draw one block waveform for one tree entry, with CFD and integration window.
// ---------------------------------------------------------------------------
void draw_hcal_waveform(const char* infile,
                        Long64_t entry,
                        int blockID,
                        double cfdFrac = 0.30,
                        int nPre = 4,
                        int nPost = 12,
                        int nPed = 8)
{
  using namespace HCalWF;

  TTree* T = nullptr;
  std::unique_ptr<TFile> fin(open_file_and_tree(infile, T));
  if (!fin || !T) return;
  if (!has_required_branches(T)) return;

  TTreeReader reader(T);
  TTreeReaderArray<Double_t> samps(reader, "sbs.hcal.samps");
  TTreeReaderArray<Double_t> sidx(reader, "sbs.hcal.samps_idx");
  TTreeReaderArray<Double_t> nsamps(reader, "sbs.hcal.nsamps");
  TTreeReaderArray<Double_t> elemID(reader, "sbs.hcal.samps_elemID");

  if (reader.SetEntry(entry) != TTreeReader::kEntryValid) {
    std::cerr << "ERROR: could not read entry " << entry << "\n";
    return;
  }

  const int iw = find_waveform_index(elemID, blockID);
  if (iw < 0) {
    std::cerr << "No waveform found for HCal block " << blockID
              << " in entry " << entry << "\n";
    return;
  }

  if (iw >= int(sidx.GetSize()) || iw >= int(nsamps.GetSize())) {
    std::cerr << "Waveform indexing arrays have inconsistent sizes.\n";
    return;
  }

  const int start = int(std::llround(sidx[iw]));
  const int n = int(std::llround(nsamps[iw]));

  PulseResult r = analyze_waveform(samps, start, n, blockID,
                                   cfdFrac, nPre, nPost, nPed, 0.0);
  if (!r.ok) {
    std::cerr << "Could not reconstruct waveform for block " << blockID << "\n";
    return;
  }

  std::vector<double> x(n), y(n);
  for (int i = 0; i < n; ++i) {
    x[i] = i;
    y[i] = samps[start+i] - r.baseline;
  }

  TCanvas* c = new TCanvas(Form("c_hcal_wf_%lld_%d", entry, blockID),
                           "HCal waveform", 1000, 700);
  TGraph* g = new TGraph(n, x.data(), y.data());
  g->SetTitle(Form("Entry %lld, HCal block %d;Sample;Calibrated ADC sample - residual baseline",
                   entry, blockID));
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.7);
  g->Draw("APL");

  const double ymin = g->GetYaxis()->GetXmin();
  const double ymax = g->GetYaxis()->GetXmax();

  TLine* lcfd = new TLine(r.tCFD, ymin, r.tCFD, ymax);
  lcfd->SetLineStyle(2);
  lcfd->SetLineWidth(2);
  lcfd->Draw();

  TLine* llo = new TLine(r.intLo, ymin, r.intLo, ymax);
  llo->SetLineStyle(3);
  llo->Draw();

  TLine* lhi = new TLine(r.intHi, ymin, r.intHi, ymax);
  lhi->SetLineStyle(3);
  lhi->Draw();

  TLegend* leg = new TLegend(0.58, 0.67, 0.88, 0.88);
  leg->AddEntry(g, Form("peak = %.3g, S/N = %.2f", r.peak, r.snr), "lp");
  leg->AddEntry(lcfd, Form("CFD %.0f%%: t = %.2f samples", 100*cfdFrac, r.tCFD), "l");
  leg->AddEntry(llo, Form("window [%d,%d], Q = %.4g", r.intLo, r.intHi, r.qFixed), "l");
  if (std::isfinite(r.rise1090))
    leg->AddEntry((TObject*)nullptr, Form("10-90 rise = %.2f samples", r.rise1090), "");
  leg->Draw();

  std::cout << "Entry " << entry << ", block " << blockID << "\n"
            << "  waveform index = " << iw << "\n"
            << "  start           = " << start << "\n"
            << "  nsamps          = " << n << "\n"
            << "  residual base   = " << r.baseline << "\n"
            << "  baseline RMS    = " << r.noiseRMS << "\n"
            << "  peak            = " << r.peak << "\n"
            << "  S/N             = " << r.snr << "\n"
            << "  CFD time        = " << r.tCFD << " samples\n"
            << "  fixed Q         = " << r.qFixed << "\n"
            << "  full Q          = " << r.qFull << "\n";
}


// ---------------------------------------------------------------------------
// Reconstruct the current PRIMARY CLUSTER using custom waveform integration.
// Writes an event-level output tree and diagnostic histograms.
// ---------------------------------------------------------------------------
void hcal_waveform_energy_study(const char* infile,
                                const char* outname = "hcal_waveform_energy_study.root",
                                Long64_t nEvents = -1,
                                double cfdFrac = 0.30,
                                int nPre = 4,
                                int nPost = 12,
                                int nPed = 8,
                                double minBlockSNR = 0.0,
                                double minReplayE = 0.0,
                                double maxReplayE = 1.0e9)
{
  using namespace HCalWF;

  TTree* T = nullptr;
  std::unique_ptr<TFile> fin(open_file_and_tree(infile, T));
  if (!fin || !T) return;
  if (!has_required_branches(T)) return;

  TTreeReader reader(T);
  TTreeReaderArray<Double_t> samps(reader, "sbs.hcal.samps");
  TTreeReaderArray<Double_t> sidx(reader, "sbs.hcal.samps_idx");
  TTreeReaderArray<Double_t> nsamps(reader, "sbs.hcal.nsamps");
  TTreeReaderArray<Double_t> elemID(reader, "sbs.hcal.samps_elemID");
  TTreeReaderArray<Double_t> clusID(reader, "sbs.hcal.clus_blk.id");
  TTreeReaderValue<Double_t> replayE(reader, "sbs.hcal.e");

  TFile* fout = new TFile(outname, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: could not create " << outname << "\n";
    return;
  }

  Long64_t o_entry = -1;
  Int_t o_nclusblk = 0;
  Int_t o_nused = 0;
  Int_t o_nmissing = 0;
  Double_t o_replayE = 0.0;
  Double_t o_qFixed = 0.0;
  Double_t o_qFull = 0.0;
  Double_t o_ratio = 0.0;
  Double_t o_meanNoiseRMS = 0.0;
  Double_t o_meanRise1090 = 0.0;
  Double_t o_meanRise2080 = 0.0;

  TTree* W = new TTree("W", "Custom HCal waveform reconstruction");
  W->Branch("entry", &o_entry, "entry/L");
  W->Branch("nclusblk", &o_nclusblk, "nclusblk/I");
  W->Branch("nused", &o_nused, "nused/I");
  W->Branch("nmissing", &o_nmissing, "nmissing/I");
  W->Branch("replayE", &o_replayE, "replayE/D");
  W->Branch("qFixed", &o_qFixed, "qFixed/D");
  W->Branch("qFull", &o_qFull, "qFull/D");
  W->Branch("ratio", &o_ratio, "ratio/D");
  W->Branch("meanNoiseRMS", &o_meanNoiseRMS, "meanNoiseRMS/D");
  W->Branch("meanRise1090", &o_meanRise1090, "meanRise1090/D");
  W->Branch("meanRise2080", &o_meanRise2080, "meanRise2080/D");

  TH1D* hRatio = new TH1D("hRatio",
      "Custom primary-cluster charge / replay HCal energy;Q_{custom}/E_{replay};Events",
      300, 0.0, 3.0);
  TH2D* hCorr = new TH2D("hCorr",
      "Custom waveform charge vs replay HCal energy;Replay E;Custom Q",
      200, 0.0, 1.0, 200, 0.0, 1.0);
  hCorr->SetCanExtend(TH1::kAllAxes);
  TH1D* hNoise = new TH1D("hNoise",
      "Residual baseline RMS;RMS of calibrated samples;Blocks", 300, 0.0, 50.0);
  TH1D* hSNR = new TH1D("hSNR", "Waveform peak S/N;peak/RMS;Blocks", 300, 0.0, 300.0);
  TH1D* hRise1090 = new TH1D("hRise1090",
      "10-90% leading-edge rise time;Rise time (samples);Blocks", 120, 0.0, 12.0);
  TH1D* hRise2080 = new TH1D("hRise2080",
      "20-80% leading-edge rise time;Rise time (samples);Blocks", 120, 0.0, 12.0);
  TH1D* hMissing = new TH1D("hMissing",
      "Primary-cluster blocks without matching stored waveform;N missing blocks/event;Events",
      10, -0.5, 9.5);

  Long64_t maxEntries = T->GetEntries();
  if (nEvents >= 0) maxEntries = std::min(maxEntries, nEvents);

  Long64_t nread = 0;
  Long64_t naccepted = 0;
  Long64_t totalMissing = 0;

  while (reader.Next()) {
    if (nread >= maxEntries) break;
    const Long64_t entry = reader.GetCurrentEntry();
    ++nread;

    const double erep = *replayE;
    if (!(erep > minReplayE && erep < maxReplayE)) continue;
    if (clusID.GetSize() == 0) continue;

    double qFixedSum = 0.0;
    double qFullSum = 0.0;
    double noiseSum = 0.0;
    double rise10Sum = 0.0;
    double rise20Sum = 0.0;
    int nNoise = 0;
    int nRise10 = 0;
    int nRise20 = 0;
    int nused = 0;
    int nmissing = 0;

    for (int ib = 0; ib < int(clusID.GetSize()); ++ib) {
      const int id = int(std::llround(clusID[ib]));
      const int iw = find_waveform_index(elemID, id);
      if (iw < 0 || iw >= int(sidx.GetSize()) || iw >= int(nsamps.GetSize())) {
        ++nmissing;
        continue;
      }

      const int start = int(std::llround(sidx[iw]));
      const int n = int(std::llround(nsamps[iw]));
      PulseResult p = analyze_waveform(samps, start, n, id,
                                       cfdFrac, nPre, nPost, nPed, minBlockSNR);
      if (!p.ok) continue;

      qFixedSum += p.qFixed;
      qFullSum += p.qFull;
      noiseSum += p.noiseRMS;
      ++nNoise;
      ++nused;

      hNoise->Fill(p.noiseRMS);
      if (std::isfinite(p.snr)) hSNR->Fill(p.snr);
      if (std::isfinite(p.rise1090)) {
        hRise1090->Fill(p.rise1090);
        rise10Sum += p.rise1090;
        ++nRise10;
      }
      if (std::isfinite(p.rise2080)) {
        hRise2080->Fill(p.rise2080);
        rise20Sum += p.rise2080;
        ++nRise20;
      }
    }

    if (nused <= 0) continue;

    o_entry = entry;
    o_nclusblk = int(clusID.GetSize());
    o_nused = nused;
    o_nmissing = nmissing;
    o_replayE = erep;
    o_qFixed = qFixedSum;
    o_qFull = qFullSum;
    o_ratio = (erep != 0.0) ? qFixedSum/erep : 0.0;
    o_meanNoiseRMS = (nNoise > 0) ? noiseSum/nNoise : 0.0;
    o_meanRise1090 = (nRise10 > 0) ? rise10Sum/nRise10
                                    : std::numeric_limits<double>::quiet_NaN();
    o_meanRise2080 = (nRise20 > 0) ? rise20Sum/nRise20
                                    : std::numeric_limits<double>::quiet_NaN();

    W->Fill();
    hRatio->Fill(o_ratio);
    hCorr->Fill(erep, qFixedSum);
    hMissing->Fill(nmissing);
    totalMissing += nmissing;
    ++naccepted;
  }

  fout->cd();
  W->Write();
  hRatio->Write();
  hCorr->Write();
  hNoise->Write();
  hSNR->Write();
  hRise1090->Write();
  hRise2080->Write();
  hMissing->Write();
  fout->Close();

  std::cout << "\nHCal waveform study finished\n"
            << "  entries read       = " << nread << "\n"
            << "  accepted events    = " << naccepted << "\n"
            << "  missing waveforms  = " << totalMissing << "\n"
            << "  CFD fraction       = " << cfdFrac << "\n"
            << "  integration window = -" << nPre << " / +" << nPost
            << " samples around floor(t_CFD)\n"
            << "  pedestal-tail bins = " << nPed << "\n"
            << "  minimum block S/N  = " << minBlockSNR << "\n"
            << "  input              = " << infile << "\n"
            << "  output             = " << outname << "\n";

  if (totalMissing > 0) {
    std::cout << "\nNOTE: Some primary-cluster block IDs did not match samps_elemID.\n"
              << "      If the missing fraction is large, check whether your replay\n"
              << "      uses a 0/1-based block-ID convention difference.\n";
  }
}


// ---------------------------------------------------------------------------
// Scan CFD fraction and integration-window sizes.
//
// Metric = fractional RMS of Qcustom / sbs.hcal.e across selected events.
// This compares custom reconstruction stability to the existing replay energy.
// It is NOT yet an independent physical energy-resolution measurement.
// ---------------------------------------------------------------------------
void hcal_scan_windows(const char* infile,
                       Long64_t nEvents = -1,
                       double minReplayE = 0.0,
                       double maxReplayE = 1.0e9,
                       int nPed = 8,
                       double minBlockSNR = 0.0)
{
  using namespace HCalWF;

  TTree* T = nullptr;
  std::unique_ptr<TFile> fin(open_file_and_tree(infile, T));
  if (!fin || !T) return;
  if (!has_required_branches(T)) return;

  const std::vector<double> cfdList = {0.20, 0.30, 0.40, 0.50};
  const std::vector<int> preList = {2, 3, 4, 5, 6};
  const std::vector<int> postList = {8, 10, 12, 14, 16, 18};

  struct Stat {
    Long64_t n = 0;
    double sum = 0.0;
    double sum2 = 0.0;
  };

  const int nCFD = int(cfdList.size());
  const int nPreN = int(preList.size());
  const int nPostN = int(postList.size());
  std::vector<Stat> stats(nCFD*nPreN*nPostN);

  auto index = [=](int ic, int ipre, int ipost) {
    return (ic*nPreN + ipre)*nPostN + ipost;
  };

  TTreeReader reader(T);
  TTreeReaderArray<Double_t> samps(reader, "sbs.hcal.samps");
  TTreeReaderArray<Double_t> sidx(reader, "sbs.hcal.samps_idx");
  TTreeReaderArray<Double_t> nsamps(reader, "sbs.hcal.nsamps");
  TTreeReaderArray<Double_t> elemID(reader, "sbs.hcal.samps_elemID");
  TTreeReaderArray<Double_t> clusID(reader, "sbs.hcal.clus_blk.id");
  TTreeReaderValue<Double_t> replayE(reader, "sbs.hcal.e");

  Long64_t maxEntries = T->GetEntries();
  if (nEvents >= 0) maxEntries = std::min(maxEntries, nEvents);

  Long64_t nread = 0;

  while (reader.Next()) {
    if (nread >= maxEntries) break;
    ++nread;

    const double erep = *replayE;
    if (!(erep > minReplayE && erep < maxReplayE)) continue;
    if (clusID.GetSize() == 0) continue;

    // Store waveform start/length/ID once for the current primary cluster.
    struct WInfo { int id; int start; int n; };
    std::vector<WInfo> winfo;
    bool missing = false;

    for (int ib = 0; ib < int(clusID.GetSize()); ++ib) {
      const int id = int(std::llround(clusID[ib]));
      const int iw = find_waveform_index(elemID, id);
      if (iw < 0 || iw >= int(sidx.GetSize()) || iw >= int(nsamps.GetSize())) {
        missing = true;
        break;
      }
      WInfo w;
      w.id = id;
      w.start = int(std::llround(sidx[iw]));
      w.n = int(std::llround(nsamps[iw]));
      winfo.push_back(w);
    }
    if (missing || winfo.empty()) continue;

    for (int ic = 0; ic < nCFD; ++ic) {
      for (int ipre = 0; ipre < nPreN; ++ipre) {
        for (int ipost = 0; ipost < nPostN; ++ipost) {

          double qsum = 0.0;
          int nused = 0;

          for (const WInfo& w : winfo) {
            PulseResult p = analyze_waveform(samps, w.start, w.n, w.id,
                                             cfdList[ic], preList[ipre], postList[ipost],
                                             nPed, minBlockSNR);
            if (!p.ok) continue;
            qsum += p.qFixed;
            ++nused;
          }

          if (nused != int(winfo.size())) continue;
          const double ratio = qsum / erep;
          if (!std::isfinite(ratio)) continue;

          Stat& s = stats[index(ic,ipre,ipost)];
          ++s.n;
          s.sum += ratio;
          s.sum2 += ratio*ratio;
        }
      }
    }
  }

  double bestFracRMS = std::numeric_limits<double>::infinity();
  double bestCFD = 0.0;
  int bestPre = 0;
  int bestPost = 0;
  Long64_t bestN = 0;
  double bestMean = 0.0;
  double bestSigma = 0.0;

  std::cout << "\nCFD/window scan results\n"
            << " metric = RMS(Qcustom/Ereplay) / mean(Qcustom/Ereplay)\n\n"
            << " CFD    pre  post      N       mean ratio     sigma      fracRMS\n"
            << " ------------------------------------------------------------------\n";

  for (int ic = 0; ic < nCFD; ++ic) {
    for (int ipre = 0; ipre < nPreN; ++ipre) {
      for (int ipost = 0; ipost < nPostN; ++ipost) {
        const Stat& s = stats[index(ic,ipre,ipost)];
        if (s.n < 2) continue;

        const double m = s.sum / double(s.n);
        const double var = std::max(0.0, s.sum2/double(s.n) - m*m);
        const double sig = std::sqrt(var);
        const double frac = (std::fabs(m) > 0.0) ? sig/std::fabs(m)
                                                 : std::numeric_limits<double>::infinity();

        std::cout << Form(" %.2f    %2d    %2d   %7lld    %10.5g   %9.4g   %9.5f\n",
                          cfdList[ic], preList[ipre], postList[ipost],
                          s.n, m, sig, frac);

        if (frac < bestFracRMS) {
          bestFracRMS = frac;
          bestCFD = cfdList[ic];
          bestPre = preList[ipre];
          bestPost = postList[ipost];
          bestN = s.n;
          bestMean = m;
          bestSigma = sig;
        }
      }
    }
  }

  if (std::isfinite(bestFracRMS)) {
    std::cout << "\nBest reconstruction-consistency setting:\n"
              << "  CFD fraction = " << bestCFD << "\n"
              << "  pre samples  = " << bestPre << "\n"
              << "  post samples = " << bestPost << "\n"
              << "  N events     = " << bestN << "\n"
              << "  mean ratio   = " << bestMean << "\n"
              << "  sigma ratio  = " << bestSigma << "\n"
              << "  frac. RMS    = " << bestFracRMS << "\n";
  } else {
    std::cout << "No valid scan points were found.\n";
  }

  std::cout << "\nIMPORTANT: This scan keeps the ORIGINAL SBS primary-cluster membership.\n"
            << "It optimizes waveform timing/integration only; it does not yet\n"
            << "rebuild clusters with new seed/neighbor thresholds.\n";
}

