bool normCB = false;
//bool normCB = true;
int nPeaks;
double binFactor;
double binF;
double rangeMin, rangeMax;
double eMin = 6200;				// lower Si alpha energy value
double eMax = 6800;				// upper Si alpha energy value
double conf_int = 0.683;

double fGam(double *x, double *par);
double   bg(double *x, double *par);
double   CB(double *x, double *par);
double CBBG(double *x, double *par);
double  CB1(double *x, double *par);
TH1F *getHisto(const char *histo);
