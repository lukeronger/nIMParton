#ifndef nIMPARTON_H_
#define nIMPARTON_H_ 1

class nIMParton
//this class is for users
{
public:
	nIMParton(unsigned int Z_temp=1, unsigned int A_temp=1); //constructor function, the first parameter is Z, the other is A for a nucleus
	virtual void setDataSet(int);                            //Choose a data set, 1 is for set A and 2 is for set B
	virtual double getRToN(int, double x, double Q2) const;  //Get parton ratio of nucleus to free nucleon
	virtual double getRToD(int, double x, double Q2) const;  //Get parton ratio of nucleus to deuteron
	virtual double getF2RToN(double x, double Q2) const;     //Get F2 ratio of nucleus to free nucleon
	virtual double getF2RToD(double x, double Q2) const;     //Get F2 ratio of nucleus to deuteron
	virtual ~nIMParton(void);                                //deconstructor function

private:
	unsigned int xMax;         //grid points number for variable x
	unsigned int Q2Max;        //grid points number for variable Q^2
	unsigned int flavorMax;    //flavor number in the grid data
	double lnxstep;            //step length of ln(x) for the grid data
	double lnQ2step;           //step length of ln(Q^2) for the grid data
	unsigned int Z;            //atomic number for a nucleus
	unsigned int A;            //mass number for a nucleus

        double * gridN;            //grid data array for free nucleon
        double * gridD;            //grid data array for deuteron
        double * gridDA;           //data array storing data set A for deuteron
        double * gridDB;           //data array storing data set B for deuteron
        double * gridA;            //grid data array for a nucleus
        double * gridAA;           //data array storing data set A for a nucleus
        double * gridAB;           //data array storing data set B for a nucleus

	virtual double getPDF(int, double x, double Q2) const;   //user function to get PDFs of a nucleus, see ReadMe.txt for details   
        virtual double getXUV(double x, double Q2) const;        //return x(u -ubaar) of a nucleus
        virtual double getXDV(double x, double Q2) const;        //return x(d - dbar) of a nucleus
        virtual double getXUSea(double x, double Q2) const;      //return 2x*ubar of a nucleus
        virtual double getXDSea(double x, double Q2) const;      //return 2x*dbar of a nucleus
        virtual double getXSSea(double x, double Q2) const;      //return 2x*sbar of a nucleus
        virtual double getF2(double x, double Q2) const;         //return structure function F2 of a nucleus
        virtual double getXGluon(double x, double Q2) const;     //return x*gluon of a nucleus
	virtual double getPDFN(int, double x, double Q2) const;  //user function to get PDFs of free nucleon, see ReadMe.txt for details       
        virtual double getXUVN(double x, double Q2) const;       //return x(u -ubaar) of free nucleon
        virtual double getXDVN(double x, double Q2) const;       //return x(d - dbar) of free nucleon
        virtual double getXUSeaN(double x, double Q2) const;     //return 2x*ubar of free nucleon
        virtual double getXDSeaN(double x, double Q2) const;     //return 2x*dbar of free nucleon
        virtual double getXSSeaN(double x, double Q2) const;     //return 2x*sbar of free nucleon
        virtual double getF2N(double x, double Q2) const;        //return structure function F2 of free nucleon 
        virtual double getXGluonN(double x, double Q2) const;    //return x*gluon of free nucleon
	virtual double getPDFD(int, double x, double Q2) const;  //user function to get PDFs of deuteron, see ReadMe.txt for details       
        virtual double getXUVD(double x, double Q2) const;       //return x(u -ubaar) of deuteron
        virtual double getXDVD(double x, double Q2) const;       //return x(d - dbar) of deuteron
        virtual double getXUSeaD(double x, double Q2) const;     //return 2x*ubar of deuteron
        virtual double getXDSeaD(double x, double Q2) const;     //return 2x*dbar of deuteron
        virtual double getXSSeaD(double x, double Q2) const;     //return 2x*sbar of deuteron
        virtual double getF2D(double x, double Q2) const;        //return structure function F2 of deuteron
        virtual double getXGluonD(double x, double Q2) const;    //return x*gluon of deuteron

	double fitLinear(double x, double* px, double* pf) const;       //linear interpolation function
	double fitQuadratic(double x, double* px, double* pf) const;    //quadratic interpolation function
        double getPDFType(int, double x, double Q2) const;
        double getPDFTypeN(int, double x, double Q2) const;
        double getPDFTypeD(int, double x, double Q2) const; 
};

#endif


