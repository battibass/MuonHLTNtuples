#ifndef CMSSWCiemat_MuTree_H
#define CMSSWCiemat_MuTree_H

#include "TROOT.h"
#include "TMath.h"
#include <vector>
#include <string>


enum part_idx { idx_top1 = 0, idx_w1, idx_b1,
    idx_wd11, idx_wd12,
    idx_top2, idx_w2, idx_b2,
    idx_wd21, idx_wd22 };

namespace ciemat {


    enum type { VETOED = 0, MUON, ELECTRON, DIMUON, DIELECTRON, MUONELE };

    class GenInfo {
        public:
            Float_t trueNumberOfInteractions; // Number of simultaneous interactions generated
            Float_t PUweight; // Weight for MC in order to match data if sigma(minBiasInel)=68 mb
            Float_t PUweightDwForSystematics; // Weight for MC in order to match data if sigma(minBiasInel)=73.5 mb
            Float_t PUweightUpForSystematics; // Weight for MC in order to match data if sigma(minBiasInel)=73.5 mb
            Int_t id1; // Id of first parton in the hard interaction (PDGId for q, 0 for gluon)
            Int_t id2; // Id of second parton in the hard interaction (PDGId for q, 0 for gluon)
            Float_t x1; // x proton fraction of first parton in the hard interaction
            Float_t x2; // x proton fraction of second parton in the hard interaction
            Float_t scalePDF; // Energy scale in the hard interaction [GeV]
            Int_t NUP; // number of particles in the LHE diagram ( = (5+ #jets) in MadGraph V+jets)

            GenInfo(){};
            virtual ~GenInfo(){};

            ClassDef(GenInfo,1)
    };


    class TTGenInfo {
        public:
            Float_t cosThetaStar;     // 
            Float_t cosThetaStar1;    //
            Float_t cosThetaStarHad;  // 
            Float_t cosThetaStarHad1; // 
            Float_t cosThetaStarHadWrong; // wrong q choice
            Bool_t isSemiLeptonic;
            Bool_t isFullLeptonic;
            Bool_t isFullHadronic;
            Bool_t isMuon;
            Bool_t isElectron;
            Bool_t isTau;
            Bool_t buggy;

            TTGenInfo() :
                cosThetaStar(-666),
                cosThetaStar1(-666),
                cosThetaStarHad(-666),
                cosThetaStarHad1(-666),
                cosThetaStarHadWrong(-666),
                isSemiLeptonic(false),
                isFullLeptonic(false),
                isFullHadronic(false),
                isMuon(false),
                isElectron(false),
                isTau(false),
                buggy(false) {};
            virtual ~TTGenInfo(){};

            ClassDef(TTGenInfo,1)
    };

    class GenParticle {
        public:
            Int_t pdgId; // PDG identifier
            Int_t status; // MC status
            Float_t energy; // energy [GeV]
            Float_t pt; // pt [GeV]
            Float_t eta; // eta
            Float_t phi; // phi
            Float_t vx; // x coordinate of production vertex [cm]
            Float_t vy;// y coordinate of production vertex [cm]
            Float_t vz;// z coordinate of production vertex [cm]
            std::vector<Int_t> mothers; // vector of indices of mothers

            GenParticle(){};
            virtual ~GenParticle(){};

            ClassDef(GenParticle,1)
    };


    class Muon {
        public:
            Float_t pt;  // pt [GeV]
            Float_t eta; // eta
            Float_t phi; // phi

            Float_t energy ;  // GC TMP
            Float_t resolEt;  // pt [GeV]
            Float_t resolEta; // eta
            Float_t resolPhi; // phi
            Float_t resolE ;
            Float_t chargedHadronIso;
            Float_t photonIso;
            Float_t neutralHadronIso;
            Int_t   isGlobal;
            Int_t   isTight;

            Int_t   charge;    // charge
            Float_t iso_pflow; // PF isolation in dR<0.4 cone dBeta
            Float_t dxy;       // signed transverse distance to primary vertex [cm]
            Float_t dz;        // signed longitudinal distance to primary vertex at min. transv. distance [cm]
            Float_t edxy;      // uncertainty on dxy [cm]
            Float_t edz;       // uncertainty on dz [cm]
            Float_t dxybs;     // signed transverse distance to beamspot [cm]
            Float_t dzbs;      // signed longitudinal distance to beamspot [cm]

            Muon(){};//{};
            virtual ~Muon(){};

            // CB Top muon selections (do not include pt cut;  Real top selection is Iso<0.12)
            bool isTopTight() const { return fabs(eta)<2.1 && isTight && fabs(iso_pflow) < 0.125; };
            bool isTopLoose() const { return fabs(eta)<2.4 && fabs(iso_pflow) < 0.20; };

            ClassDef(Muon,1)
    };

    class Electron {
        public:
            Float_t pt; // pt [GeV]
            Float_t eta; // eta
            Float_t phi; // phi
            Float_t energy ; //
            Float_t resolEt; // pt [GeV]
            Float_t resolEta; // eta
            Float_t resolPhi; // phi
            Float_t resolE ;
            Float_t chargedHadronIso;
            Float_t photonIso;
            Float_t neutralHadronIso;
            Bool_t isNotCrackEBEE; // !isEBEEGap
            Bool_t isNotCrack; // 1.4442 > |superCluster.eta| > 1.5560
            Bool_t passConvRej;

            // Lepton Id
            Bool_t isVetoForDil;
            Bool_t isVetoForLj;
            Bool_t isTightForDil;
            Bool_t isTightForLj;

            Float_t mvaTrigV0;
            Float_t mvaNonTrigV0;
            Bool_t cbTight;
            Bool_t cbVeto;
            // END GC TMP

            Int_t charge; // charge
            Float_t iso_pflow; // PF isolation in dR<0.3 cone, EA-subtracted
            Float_t iso_pflow_base; // PF isolation in dR<0.3 cone
            Float_t dxy; // signed transverse distance to primary vertex [cm]
            Float_t dz; // signed longitudinal distance to primary vertex at min. transv. distance [cm]
            Float_t edxy; // uncertainty on dxy [cm]
            Float_t edz; // uncertainty on dz [cm]
            Float_t dxybs; // signed transverse distance to beamspot [cm]
            Float_t dzbs; // signed longitudinal distance to beamspot [cm]

            Electron(){};//{};
            virtual ~Electron(){};

            ClassDef(Electron,1)
    };

    class SSV {
        public:
            Int_t nTracks; // #tracks at this secondary vertex
            Float_t flightDistance; // Flight distance [cm] (w.r.t. PV)
            Float_t errorFlightDistance; // Errors on flight distance [cm]
            Float_t sig2d; // Transverse significance
            Float_t vertexPhi; // Phi of 3d-vector SSV-PV
            Float_t vertexEta; // Eta of 3d-vector SSV-PV

            SSV(){};
            virtual ~SSV(){};

            ClassDef(SSV,1)
    };

    class DMass {
        public:
            Float_t MassSSV; // Mass assumming all pions [GeV]
            Float_t MassD; // Mass closer to D0 or D+- assuming one Kaon [GeV]
            Float_t MassDs; // Mass closer to Ds+- assuming two Kaons [GeV] (only for odd number of tracks)

            DMass(){};
            virtual ~DMass(){};

            ClassDef(DMass,1)
    };

    class Jet {
        public:
            Float_t et; // et [GeV]
            Float_t pt; // pt [GeV]
            Float_t eta; // eta
            Float_t phi; // phi

            Float_t ptGen; // pt [GeV]


            //// GC TMP
            Float_t energy; // GC TMP
            Float_t resolEt; // pt [GeV]
            Float_t resolEta; // eta
            Float_t resolPhi; // phi
            Float_t resolE;
            Float_t resolEtB; // pt [GeV]
            Float_t resolEtaB; // eta
            Float_t resolPhiB; // phi
            Float_t resolEB;
            Int_t partonFlavour;
            // END GC TMP

            Int_t chargedMultiplicity; // charged multiplicity
            Float_t emEnergy; // electromagnetic energy [GeV]
            Float_t muonEnergy; // muonic energy [GeV]
            Float_t jetunc; // estimated jet energy scale uncertainty
            std::vector<ciemat::Muon> muons; // vector of muons in jets
            Float_t btagTCHE; // TCHE discriminant
            Float_t negativeBtagTCHE; // negative TCHE discriminant
            Float_t btagCSV; // CSV discriminant
            std::vector<ciemat::SSV> secondaryVertices; // vector of Simple Sec. Vertices (SSV)
            std::vector<ciemat::DMass> secondaryMasses; // vector of masses with sec. tracks

            Jet(){};//{};
            virtual ~Jet(){};

            double massSSV() const { // invariant mass of SSV tracks, assuming pion mass [GeV]
                if (secondaryMasses.size()>0) return secondaryMasses[0].MassSSV;
                else return -1.;
            }

            double ptmu() const { // pt of most energetic muon in jet [GeV]
                double ptmax = 0.;
                for (unsigned int im=0; im<muons.size(); ++im) {
                    double pt = muons[im].pt;
                    if (pt>ptmax) ptmax = pt;
                }
                return ptmax;
            }

            double ptrel() const { // ptrel of most energetic muon in jet [GeV]
                int index = -1;
                double ptmax = 0.;
                for (unsigned int im=0; im<muons.size(); ++im) {
                    double pt = muons[im].pt;
                    if (pt>ptmax) {
                        ptmax = pt;
                        index = im;
                    }
                }
                if (index>=0) return ptrel(index); else return 0.;
            }

            double ptrel(unsigned int im) const { // ptrel w.r.t to jet direction of muon im [GeV]
                if (im>=muons.size()) return 0.;
                ciemat::Muon mu = muons[im];
                double pxm = mu.pt*cos(mu.phi);
                double pym = mu.pt*sin(mu.phi);
                double pzm = mu.pt*sinh(mu.eta);
                double p2m = pxm*pxm+pym*pym+pzm*pzm;
                double pxj = pt*cos(phi);
                double pyj = pt*sin(phi);
                double pzj = pt*sinh(eta);
                double p2j = pxj*pxj+pyj*pyj+pzj*pzj;
                double projm = pxm*pxj+pym*pyj+pzm*pzj;
                double abs_ptrel = sqrt(p2m - projm*projm/p2j);
                if (mu.charge>0) return abs_ptrel; else return -abs_ptrel;
            }

            double numberSSV() const { // number of positive simple secondary vertices (SSVs)
                int cont=0;
                unsigned int nSSV = secondaryVertices.size();
                for (unsigned int isv=0; isv<nSSV; ++isv) {
                    const ciemat::SSV& svtx = secondaryVertices[isv];
                    if (svtx.nTracks<2) continue;
                    if (svtx.errorFlightDistance<=0.) continue;
                    if (svtx.flightDistance<0.) continue;
                    cont++;
                }
                return cont;
            }

            double numbernegativeSSV() const { // number of negative SSVs
                int cont=0;
                unsigned int nSSV = secondaryVertices.size();
                for (unsigned int isv=0; isv<nSSV; ++isv) {
                    const ciemat::SSV& svtx = secondaryVertices[isv];
                    if (svtx.nTracks<2) continue;
                    if (svtx.errorFlightDistance<=0.) continue;
                    if (svtx.flightDistance>0.) continue;
                    cont++;
                }
                return cont;
            }

            double btagSSVHE() const { // SSVHE positive discriminant
                unsigned int nSSV = secondaryVertices.size();
                for (unsigned int isv=0; isv<nSSV; ++isv) {
                    const ciemat::SSV& svtx = secondaryVertices[isv];
                    if (svtx.nTracks<2) continue;
                    if (svtx.errorFlightDistance<=0.) continue;
                    if (svtx.flightDistance<0.) continue;
                    return log(1.+svtx.flightDistance/svtx.errorFlightDistance);
                }
                return -9999.;
            };
            double btagSSVHP() const { // SSVHP positive discriminant
                unsigned int nSSV = secondaryVertices.size();
                for (unsigned int isv=0; isv<nSSV; ++isv) {
                    const ciemat::SSV& svtx = secondaryVertices[isv];
                    if (svtx.nTracks<3) continue;
                    if (svtx.errorFlightDistance<=0.) continue;
                    if (svtx.flightDistance<0.) continue;
                    return log(1.+svtx.flightDistance/svtx.errorFlightDistance);
                }
                return -9999.;
            };
            double btagNegativeSSVHE() const { // SSVHE negative discriminant
                unsigned int nSSV = secondaryVertices.size();
                for (unsigned int isv=0; isv<nSSV; ++isv) {
                    const ciemat::SSV& svtx = secondaryVertices[isv];
                    if (svtx.nTracks<2) continue;
                    if (svtx.errorFlightDistance<=0.) continue;
                    if (svtx.flightDistance>0.) continue;
                    return -log(1.-svtx.flightDistance/svtx.errorFlightDistance);
                }
                return -9999.;
            };
            double btagNegativeSSVHP() const { // SSVHP negative discriminant
                unsigned int nSSV = secondaryVertices.size();
                for (unsigned int isv=0; isv<nSSV; ++isv) {
                    const ciemat::SSV& svtx = secondaryVertices[isv];
                    if (svtx.nTracks<3) continue;
                    if (svtx.errorFlightDistance<=0.) continue;
                    if (svtx.flightDistance>0.) continue;
                    return -log(1.-svtx.flightDistance/svtx.errorFlightDistance);
                }
                return -9999.;
            };

            ClassDef(Jet,1)
    };



    //  bool passed =  event.hlt.find( "HLT_Mu15_and_whatever" );

    class HLT {
        public:
            std::vector<std::string> triggers; // vector of strings with HLT paths fired
            Float_t muMaxPt;
            Float_t eleMaxPt;

            HLT(){};
            virtual ~HLT(){};
            bool match( const std::string & path ) {
                if (  std::find (  triggers.begin(), triggers.end(), path ) != triggers.end() )
                    return true;

                return false;
            }
            bool find( const std::string & path ) {
                for ( std::vector<std::string>::const_iterator it = triggers.begin(); it != triggers.end(); ++it ) {
                    if ( it->find ( path ) != std::string::npos ) return true;
                }
                return false;
            }
            ClassDef(HLT,1)
    };



    class MET {
        public:
            Float_t met;          // MET [GeV]
            Float_t phi;          // phi of MET
            Float_t sumEt;        // Sum of transverse energies in the event [GeV]
            Float_t significance; // MET significance

            UInt_t filterResults; // word for bitmap of filter results

            Float_t resolEt;      // pt [GeV]
            Float_t resolPhi;     // phi

            MET(){};
            virtual ~MET(){};

            ClassDef(MET,1)
    };

    class Event {
        public:
            Int_t runNumber; // run number
            Int_t luminosityBlockNumber; // luminosity block number
            Int_t eventNumber; // event number
            Int_t nvvertex; // number of valid reconstructed primary vertices 
            Float_t Rho; // dEnergy/dR density from pileup+underlying activity [GeV]
            Float_t primaryVertex[3]; // 3d coordinates of PV [cm]
            Float_t cov_primaryVertex[3][3]; // 3x3 covariance matrix of PV estimation [cm*cm]
            std::vector <ciemat::GenInfo> genInfos; // venctor of genInfosi; size=0 in data
            std::vector<ciemat::GenParticle> genParticles; // venctor of genParticles size=0 in data
            std::vector<ciemat::Muon> muons; // vector of muons
            std::vector<ciemat::Electron> electrons; // vector of electrons
            std::vector<ciemat::Jet> jets; // vector of jets
            ciemat::HLT hlt; // HLT object
            ciemat::MET met; // MET object
            std::vector<ciemat::TTGenInfo> ttGenInfo; // vector of jets

            Event(){};
            virtual ~Event(){};

            GenParticle & GenTop1()   { return genParticles[idx_top1] ; };
            GenParticle & GenW1()     { return genParticles[idx_w1]; };
            GenParticle & GenB1()     { return genParticles[idx_b1]; };
            GenParticle & GenWD11()   { return genParticles[idx_wd11]; };
            GenParticle & GenWD12()   { return genParticles[idx_wd12]; };

            GenParticle & GenTop2()   { return genParticles[idx_top2] ; };
            GenParticle & GenW2()     { return genParticles[idx_w2]; };
            GenParticle & GenB2()     { return genParticles[idx_b2]; };
            GenParticle & GenWD21()   { return genParticles[idx_wd21]; };
            GenParticle & GenWD22()   { return genParticles[idx_wd22]; };


            Float_t mt(const Muon& mu) const { // MT from muon mu + MET
                return sqrt(2*mu.pt*met.met*(1-cos(mu.phi-met.phi)));
            }
            Float_t mt(const Electron& el) const { // MT from electron el + MET
                return sqrt(2*el.pt*met.met*(1-cos(el.phi-met.phi)));
            }
            Float_t acop(const Muon& mu) const { // acoplanarity between muon mu and MET
                return acos(-cos(mu.phi-met.phi));
            }
            Float_t acop(const Electron& el) const { // acoplanarity between electron el and MET
                return acos(-cos(el.phi-met.phi));
            }
            Float_t mtLeadingMuon() const { // MT from leading muon + MET
                return muons.size()>0 ? mt(muons[0]) : 0.;
            }
            Float_t mtLeadingElectron() const { // MT from leading electron + MET
                return electrons.size()>0 ? mt(electrons[0]) : 0.;
            }





            ///******************* VETO **********************///
            /////// TO BE REVIEWED
            type analysisType() {

                size_t nEleVetoForDil = 0;
                size_t nEleVetoForLj = 0;
                size_t nEleTightForDil = 0;
                size_t nEleTightForLj = 0;

                size_t nMuVeto  = 0;
                size_t nMuLoose = 0;
                size_t nMuTight = 0;

                std::vector<ciemat::Electron>::const_iterator elend = electrons.end();
                for ( std::vector<ciemat::Electron>::const_iterator
                        el = electrons.begin();  el != elend; ++el ) {
                    if ( el->isVetoForDil )  ++nEleVetoForDil;
                    if ( el->isVetoForLj )   ++nEleVetoForLj;
                    if ( el->isTightForDil ) ++nEleTightForDil;
                    if ( el->isTightForLj )  ++nEleTightForLj;
                }


                std::vector<ciemat::Muon>::const_iterator muend = muons.end();
                for ( std::vector<ciemat::Muon>::const_iterator
                        mu = muons.begin();  mu != muend; ++mu ) {
                    if ( mu->isTopTight() ) ++nMuTight;
                    if ( mu->isTopLoose() ) ++nMuLoose;
                    ++nMuVeto;
                }


                if ( nEleTightForLj == 1 && nEleVetoForLj == 1 && !nMuVeto ) return ELECTRON;
                else if ( nEleTightForDil == 2 && nEleVetoForDil == 2 && !nMuVeto ) return DIELECTRON;
                else if ( nMuTight == 1 && nMuVeto == 1 && !nEleVetoForLj ) return MUON;
                else if ( nMuLoose == 2 && nMuVeto == 2 && !nEleVetoForDil ) return DIMUON;
                else if ( nMuLoose == 1 && nMuVeto == 1 && nEleTightForDil == 1 && nEleVetoForDil == 1 ) return MUONELE;
                else return VETOED;

            }


            ClassDef(Event,1)
    };

    class Summary {
        public:
            Int_t nMuons;
            Int_t nElectrons;
            Int_t nJets;
            Int_t nJetsPt30;
            Float_t met;
            Float_t bestCSV;
            Float_t secondCSV;
            Float_t maxMT;

            Summary(){};
            virtual ~Summary(){};
            ClassDef(Summary, 1)
    };
}
#endif
