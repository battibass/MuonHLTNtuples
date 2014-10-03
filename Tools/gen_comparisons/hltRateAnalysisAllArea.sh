#! /bin/bash

VERSION="v11"

# WToMuNu
./hltRateAnalysis /data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_${VERSION}_53X_HLT701_8TeV_25PU_25ns_WToMuNu.root      -1000 
./hltRateAnalysis /data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_${VERSION}_53X_HLT701_13TeV_25PU_25ns_WToMuNu.root     -1000
./hltRateAnalysis /data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_${VERSION}_62X_HLT701_13TeV_20PU_25ns_WToMuNu.root     -1000

# QCDMu2030
./hltRateAnalysis /data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_${VERSION}_53X_HLT701_8TeV_25PU_25ns_QCDMu2030.root    -1000
./hltRateAnalysis /data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_${VERSION}_53X_HLT701_13TeV_25PU_25ns_QCDMu2030.root   -1000
./hltRateAnalysis /data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_${VERSION}_62X_HLT701_13TeV_20PU_25ns_QCDMu2030.root   -1000

# QCDMu3050
./hltRateAnalysis /data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_${VERSION}_53X_HLT701_8TeV_25PU_25ns_QCDMu3050.root    -1000
./hltRateAnalysis /data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_${VERSION}_53X_HLT701_13TeV_25PU_25ns_QCDMu3050.root   -1000
./hltRateAnalysis /data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_${VERSION}_62X_HLT701_13TeV_20PU_25ns_QCDMu3050.root   -1000

# QCDMu5080
./hltRateAnalysis /data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_${VERSION}_53X_HLT701_8TeV_25PU_25ns_QCDMu5080.root    -1000
./hltRateAnalysis /data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_${VERSION}_53X_HLT701_13TeV_25PU_25ns_QCDMu5080.root   -1000
./hltRateAnalysis /data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_${VERSION}_62X_HLT701_13TeV_20PU_25ns_QCDMu5080.root   -1000

# QCDMu80120
./hltRateAnalysis /data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_${VERSION}_53X_HLT701_8TeV_25PU_25ns_QCDMu80120.root   -1000
./hltRateAnalysis /data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_${VERSION}_53X_HLT701_13TeV_25PU_25ns_QCDMu80120.root  -1000
./hltRateAnalysis /data1/battilan/MuonHLT/RobertoRateChecks/MuonHltTree_${VERSION}_62X_HLT701_13TeV_20PU_25ns_QCDMu80120.root  -1000


