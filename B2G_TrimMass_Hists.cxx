#include "UHH2/TstarTstar/include/B2G_TrimMass_Hists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <string>

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

float inv_mass_B2G(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }


B2G_TrimMass_Hists::B2G_TrimMass_Hists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  // booking all histograms
  book<TH1F>("oneBin", "oneBin", 1, 0, 1);

  // basic variables
  book<TH1F>("mjj", "mjj", 50, 500, 2500);
  book<TH1F>("mj", "mj", 50, 0, 500);
  book<TH1F>("HT", "HT", 50, 0, 4000);

  // now with various additional requirements
  book<TH1F>("HT_abovemj_55", "HT_abovemj_55", 50, 0, 4000);
  book<TH1F>("mjj_abovemj_55", "mjj_abovemj_55", 50, 0, 3000);

  book<TH1F>("HT_abovemj_35", "HT_abovemj_35", 50, 0, 4000);
  book<TH1F>("mjj_abovemj_35", "mjj_abovemj_35", 50, 0, 3000);

  book<TH1F>("mjj_aboveHT_1000", "mjj_aboveHT_1000", 50, 0, 3000);
  book<TH1F>("mj_above_600pt", "mj_above_600pt", 50, 0, 500);

  book<TH1F>("mjj_abovemj_both55", "mjj_abovemj_both55", 50, 0, 3000);
  book<TH1F>("mj_abovemj_both55", "mj_abovemj_both55", 50, 0, 500);
  book<TH1F>("HT_abovemj_both55", "HT_abovemj_both55", 50, 0, 3000);

  book<TH1F>("mjj_abovemj_both55_aboveHT_1000", "mjj_abovemj_both55_aboveHT_1000", 50, 0, 3000);

}


void B2G_TrimMass_Hists::fill(const Event & event){
  // fill the histograms

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  hist("oneBin")->Fill(0.5, weight);

  double HT = 0.;
  for (const auto & jet : *event.topjets) HT += jet.pt();

  double mj = 0;
  if(event.topjets->size() > 0) mj = event.topjets->at(0).softdropmass();

  double mjj = 0;
  if(event.topjets->size() > 1) mjj = inv_mass_B2G(event.topjets->at(0).v4() + event.topjets->at(1).v4());

  // mj hists
  if(event.topjets->size() > 0) {
    hist("mj")->Fill(mj, weight);
    if(event.topjets->at(0).pt() > 600) hist("mj_above_600pt")->Fill(mj, weight);
  }

  // mjj hists
  if(event.topjets->size() > 1) {
    hist("mjj")->Fill(mjj, weight);
    if(mj > 35) hist("mjj_abovemj_35")->Fill(mjj, weight);
    if(mj > 55) hist("mjj_abovemj_55")->Fill(mjj, weight);
    if(HT > 1000) hist("mjj_aboveHT_1000")->Fill(mjj, weight);
  }

  // HT hists
  if(event.topjets->size() > 0) {
    hist("HT")->Fill(HT, weight);
    if(mj > 55) hist("HT_abovemj_55")->Fill(HT, weight);
    if(mj > 35) hist("HT_abovemj_35")->Fill(HT, weight);
  }

  if(event.topjets->size() > 1) {
    if((event.topjets->at(0).softdropmass() > 55) && (event.topjets->at(1).softdropmass() > 55)) {
      hist("mjj_abovemj_both55")->Fill(mjj, weight);
      hist("HT_abovemj_both55")->Fill(HT, weight);
      hist("mj_abovemj_both55")->Fill(mj, weight);
      if(HT > 1000) hist("mjj_abovemj_both55_aboveHT_1000")->Fill(mjj, weight);
    }
  }

}

B2G_TrimMass_Hists::~B2G_TrimMass_Hists(){}
