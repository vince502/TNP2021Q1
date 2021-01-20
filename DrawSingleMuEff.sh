#!/bin/bash
#void SingleMuTrigEff(int nevt=-1, bool isMC = false, int kTrigSel = 14, int hiHFBinEdge = 0, int PDtype = 1, double ptcut = 3.5)
root -l -b -q 'SingleMuTrigEff.C(-1, false, 14, 0, 1)'

root -l -b -q 'SingleMuTrigEff.C(-1, false, 13, 0, 1)'
