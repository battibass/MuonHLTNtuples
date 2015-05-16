#ifndef MuonHLT_Tools_MuonHltTree_H
#define MuonHLT_Tools_MuonHltTree_H

#include "TROOT.h"
#include "TMath.h"
#include <vector>
#include <string>

namespace muon_hlt {


  class GenInfo {
  public:
    Float_t trueNumberOfInteractions;   // Number of simultaneous interactions generated (before poissonian ev by ev smearing)
    Int_t   actualNumberOfInteractions; // Number of simultaneous interactions generated (before poissonian ev by ev smearing)
    
    GenInfo(){};
    virtual ~GenInfo(){};
    
    ClassDef(GenInfo,1)
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

    Int_t   charge;    // charge

    Int_t   isGlobal;
    Int_t   isTracker;
    Int_t   isStandAlone;

    Int_t   isSoft;
    Int_t   isLoose;
    Int_t   isTight;
    Int_t   isHighPt;
    
    Float_t chargedHadronIso;
    Float_t photonIso;
    Float_t neutralHadronIso;


    Float_t isoPflow04; // PF isolation in dR<0.4 cone dBeta
    Float_t isoPflow03; // PF isolation in dR<0.3 cone dBeta

    Float_t dxy;       // signed transverse distance to primary vertex [cm]
    Float_t dz;        // signed longitudinal distance to primary vertex at min. transv. distance [cm]
    Float_t edxy;      // uncertainty on dxy [cm]
    Float_t edz;       // uncertainty on dz [cm]
    Float_t dxybs;     // signed transverse distance to beamspot [cm]
    Float_t dzbs;      // signed longitudinal distance to beamspot [cm]

    Int_t   nHitsPixel;
    Int_t   nHitsGlobal;
    Int_t   nHitsTracker;
    Int_t   nHitsStandAlone;


    Muon(){};
    virtual ~Muon(){};

    bool isTightIso() const { return fabs(eta)<2.1 && isTight && fabs(isoPflow04) < 0.125; };
    bool isLooseIso() const { return fabs(eta)<2.1 && isLoose && fabs(isoPflow04) < 0.2;   };

    ClassDef(Muon,1)
  };

  class HLTObject {
  public:

    std::string filterTag; // name of filter passed by the object
    Float_t pt;            // pt of the object passing the filter [GeV]
    Float_t eta;           // eta of the object passing the filter
    Float_t phi;           // phi of the object passing the filter
    
    HLTObject(){};
    virtual ~HLTObject(){};

    ClassDef(HLTObject,1)

  };

  class HLT {
  public:
    std::vector<std::string> triggers; // vector of strings with HLT paths
    std::vector<muon_hlt::HLTObject>   objects;  // vector of hlt objects assing filters

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

  class Event {
  public:

    Int_t runNumber;             // run number
    Int_t luminosityBlockNumber; // luminosity block number
    Int_t eventNumber;           // event number

    Int_t nVtx;                      // number of valid reconstructed primary vertices 
    Float_t primaryVertex[3];        // 3d coordinates of PV [cm]
    Float_t cov_primaryVertex[3][3]; // 3x3 covariance matrix of PV estimation [cm*cm]

    std::vector <muon_hlt::GenInfo> genInfos;        // venctor of genInfos; size=0 in data
    std::vector<muon_hlt::GenParticle> genParticles; // venctor of genParticles size=0 in data
    std::vector<muon_hlt::Muon> muons; // vector of muons
    muon_hlt::HLT hlt;                 // HLT objects

    Event(){};
    virtual ~Event(){};

    ClassDef(Event,1)
  };

}
#endif
