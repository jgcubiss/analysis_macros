// Function declarations
//void readConfigFile();
void outputConfiguration(int optionAlphaGamma, int optionGammaGamma, int optionBetaGammaGamma);

// Parameters to be read from config file
//int optionAlphaGamma=0, optionGammaGamma=0, optionBetaGammaGamma=0;
double promptWindowAlphaGamma, promptWindowGammaGamma;
double betaEnergyMin, betaEnergyMax;

// Histos that may be used
TH2D *ag_mat, *gg_mat, *bgg_mat;

int progress=0;
int crystals=16, detectors=7;

// Variables to match types used for values stored in
// ids tree in .root files
Int_t MULT;                             //Multiplicity
ULong64_t TIME_REF;                     //Time vs ref
ULong64_t TIMESTAMP;            //Timestamp
double E_Clov[16];                       //Energy HPGe
Int_t T_Clov[16];                       //Time HPGe
Int_t M_Clov;                           //Multiplicity HPGe
double E_Si[8];          //Energy SPEDE
Int_t T_Si[8];                  //Time SPEDE
Int_t M_Si;                             //Multiplicity SPEDE
double E_Beta[2];                        //Energy beta
Int_t T_Beta[2];                        //Time beta
Int_t M_Beta;                           //Multiplicity beta
Int_t E_Proton_T1[1];           //Energy PP
Int_t T_Proton_T1[1];           //Time PP
Int_t M_Proton_T1;                      //Multiplicity PP

