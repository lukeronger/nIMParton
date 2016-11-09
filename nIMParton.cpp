#include<iostream>
#include<fstream>
#include<string>
extern "C"{
#include<math.h>
}
#include "nIMParton.h"
using namespace std;


//a method used to choose a data set
void nIMParton::setDataSet(int dataset)
{
        if(dataset==1)
        {
                gridD = gridDA;
		gridA = gridAA;
                cout<<"    Using data set A."<<endl;
        }
        else if(dataset==2)
        {
                gridD = gridDB;
                gridA = gridAB;
                cout<<"    Using data set B."<<endl;
        }
        else
        {
                cout<<"!!->Error: Unknown data set."<<endl;
                cout<<"!!->Data set should be 1 or 2, for set A and B respectviely."<<endl;
        }
}

//get parton distribution ratio of a nucleus to free nucleon
double nIMParton::getRToN(int Iparton, double x, double Q2) const
{
	return getPDF(Iparton, x, Q2)/getPDFN(Iparton, x, Q2);
}

//get parton distribution ratio of a nucleus to deuteron
double nIMParton::getRToD(int Iparton, double x, double Q2) const
{
	return getPDF(Iparton, x, Q2)/getPDFD(Iparton, x, Q2);
}

//get F2 ratio of a nucleus to free nucleon
double nIMParton::getF2RToN(double x, double Q2) const
{
	return getF2(x, Q2)/getF2N(x, Q2);
}
//get F2 ratio of a nucleus to deuteron
double nIMParton::getF2RToD(double x, double Q2) const
{
	return getF2(x, Q2)/getF2D(x, Q2);
}

//the constructor and initialization
nIMParton::nIMParton(unsigned int Z_temp, unsigned int A_temp)
:Z(Z_temp),A(A_temp)            //Z and A are parameters for a nuclei
{
	cout<<"    nIMParton version - v1.0"<<endl;
	char filename[50];
	ifstream datain;
	double x, Q2;
	unsigned int i, j;
	xMax=61;
	Q2Max=32;
	flavorMax=7;
	lnxstep=log(10)/((xMax-1)/6);
	lnQ2step=log(2.0);
	gridDA = new double[Q2Max*xMax*flavorMax];
	gridDB = new double[Q2Max*xMax*flavorMax];
	gridAA = new double[Q2Max*xMax*flavorMax];
	gridAB = new double[Q2Max*xMax*flavorMax];
	gridN = new double[Q2Max*xMax*flavorMax];
	//read grid data for interpolation
	//reading data set A
	sprintf(filename,"grid_data/gridn_%d_%d_SetA.dat",Z,A);
	cout<<"    Loading "<<filename<<endl;
	datain.open(filename);
	if(!datain.good())cout<<"!!->Error: Error while opening "<<filename<<"!\n!!->grid data file not exist?"<<endl;
	else
	for(i=0;i<Q2Max;i++)
	{
		for(j=0;j<xMax;j++)
		{
			datain>>Q2>>x>>(*(gridAA+(xMax*i+j)*7))>>(*(gridAA+(xMax*i+j)*7+1))>>(*(gridAA+(xMax*i+j)*7+2))>>(*(gridAA+(xMax*i+j)*7+3))>>(*(gridAA+(xMax*i+j)*7+4))>>(*(gridAA+(xMax*i+j)*7+5))>>(*(gridAA+(xMax*i+j)*7+6));
		}
	}
	datain.close();
	sprintf(filename,"grid_data/gridn_%d_%d_SetA.dat",1,2);
	cout<<"    Loading "<<filename<<endl;
	datain.open(filename);
	if(!datain.good())cout<<"!!->Error: Error while opening "<<filename<<"!\n!!->grid data file not exist?"<<endl;
	else
	for(i=0;i<Q2Max;i++)
	{
		for(j=0;j<xMax;j++)
		{
			datain>>Q2>>x>>(*(gridDA+(xMax*i+j)*7))>>(*(gridDA+(xMax*i+j)*7+1))>>(*(gridDA+(xMax*i+j)*7+2))>>(*(gridDA+(xMax*i+j)*7+3))>>(*(gridDA+(xMax*i+j)*7+4))>>(*(gridDA+(xMax*i+j)*7+5))>>(*(gridDA+(xMax*i+j)*7+6));
		}
	}
	datain.close();
	//reading data set B
        sprintf(filename,"grid_data/gridn_%d_%d_SetB.dat",Z,A);
        cout<<"    Loading "<<filename<<endl;
        datain.open(filename);
        if(!datain.good())cout<<"!!->Error: Error while opening "<<filename<<"!\n!!->grid data file not exist?"<<endl;
        else
        for(i=0;i<Q2Max;i++)
        {
                for(j=0;j<xMax;j++)
                {
                        datain>>Q2>>x>>(*(gridAB+(xMax*i+j)*7+0))>>(*(gridAB+(xMax*i+j)*7+1))>>(*(gridAB+(xMax*i+j)*7+2))>>(*(gridAB+(xMax*i+j)*7+3))>>(*(gridAB+(xMax*i+j)*7+4))>>(*(gridAB+(xMax*i+j)*7+5))>>(*(gridAB+(xMax*i+j)*7+6));
                }
        }
        datain.close();
        sprintf(filename,"grid_data/gridn_%d_%d_SetB.dat",1,2);
        cout<<"    Loading "<<filename<<endl;
        datain.open(filename);
        if(!datain.good())cout<<"!!->Error: Error while opening "<<filename<<"!\n!!->grid data file not exist?"<<endl;
        else
        for(i=0;i<Q2Max;i++)
        {
                for(j=0;j<xMax;j++)
                {
                        datain>>Q2>>x>>(*(gridDB+(xMax*i+j)*7+0))>>(*(gridDB+(xMax*i+j)*7+1))>>(*(gridDB+(xMax*i+j)*7+2))>>(*(gridDB+(xMax*i+j)*7+3))>>(*(gridDB+(xMax*i+j)*7+4))>>(*(gridDB+(xMax*i+j)*7+5))>>(*(gridDB+(xMax*i+j)*7+6));
                }
        }
        datain.close();
        sprintf(filename,"grid_data/gridn_%d_%d.dat",1,1);
        cout<<"    Loading "<<filename<<endl;
        datain.open(filename);
        if(!datain.good())cout<<"!!->Error: Error while opening "<<filename<<"!\n!!->grid data file not exist?"<<endl;
        else
        for(i=0;i<Q2Max;i++)
        {
                for(j=0;j<xMax;j++)
                {
                        datain>>Q2>>x>>(*(gridN+(xMax*i+j)*7+0))>>(*(gridN+(xMax*i+j)*7+1))>>(*(gridN+(xMax*i+j)*7+2))>>(*(gridN+(xMax*i+j)*7+3))>>(*(gridN+(xMax*i+j)*7+4))>>(*(gridN+(xMax*i+j)*7+5))>>(*(gridN+(xMax*i+j)*7+6));
                }
        }
        datain.close();

	//the default is set B
	gridD = gridDB;
	gridA = gridAB;

}

//the deconstructor
nIMParton::~nIMParton(void)
{
	delete[] gridN;
	delete[] gridDA;
        delete[] gridDB;
	delete[] gridAA;
        delete[] gridAB;
}

//return the parton distributions of a nucleus of different kinds at x and Q^2
double nIMParton::getPDF(int Iparton, double x, double Q2) const
{
        if(Iparton==-3 || Iparton==3)return getXSSea(x,Q2)/2.0/x;
        else if(Iparton==-2)return getXDSea(x,Q2)/2.0/x;
        else if(Iparton==2)return getXDSea(x,Q2)/2.0/x+getXDV(x,Q2)/x;
        else if(Iparton==-1)return getXUSea(x,Q2)/2.0/x;
        else if(Iparton==1)return getXUSea(x,Q2)/2.0/x+getXUV(x,Q2)/x;
        else if(Iparton==0)return getXGluon(x,Q2)/x;
        else
        {
                cout<<"!!->Error: Unknown Iparton type."<<" (Iparton = "<<Iparton<<"?)"<<endl;
                cout<<"!!->Iparton should be one of these: [-3,-2, ... 2, 3]"<<endl;
                return 0;
        }

}

//a method which returns xuv of a nucleus
double nIMParton::getXUV(double x, double Q2) const
{
	return getPDFType(1,x,Q2); 
}

//a method which returns xdv of a nucleus
double nIMParton::getXDV(double x, double Q2) const
{
        return getPDFType(2,x,Q2);
}

//a method which returns xusea of a nucleus
double nIMParton::getXUSea(double x, double Q2) const
{
        return getPDFType(3,x,Q2);
}

//a method which returns xdsea of a nucleus
double nIMParton::getXDSea(double x, double Q2) const
{
        return getPDFType(4,x,Q2);
}

//a method which returns xssea of a nucleus
double nIMParton::getXSSea(double x, double Q2) const
{
        return getPDFType(5,x,Q2);
}

//a method which returns structure function F2 of a nucleus
double nIMParton::getF2(double x, double Q2) const
{
        return getPDFType(6,x,Q2);
}

//a method which returns xgluon of a nucleus
double nIMParton::getXGluon(double x, double Q2) const
{
        return getPDFType(0,x,Q2);
}

//a method which returns different types of distributions of a nucleus
double nIMParton::getPDFType(int Iparton, double x, double Q2) const
{
	if(Iparton<0 || Iparton>6)
	{
		cout<<"!!->Error: Wrong Iparton input for getPDFType(int Iparton, double x, double Q2)."<<endl;
		return 0;
	}
	else
	{
		double lnx, lnQ2;
		int i=(int)(lnx=log(x*1e6)/lnxstep);
		int j=(int)(lnQ2=log(Q2*8)/lnQ2step);
		double g0[3], g1[3], g2[3], g[3]={0};
		if(i<0)i=0;
		if(i>(int)(xMax-3))i=xMax-3;
		if(j<0)j=0;
		if(j>29)j=29;
		//avoid log(1-x) calculation in below algorithm
		if(x>0.9999)return 0.0;
		//if x>0.5, we use A(1-x)^B to do the interpolation
		else if(x>0.5)
		{
			double vln1_x[2]={log(1-exp(log(1e-6)+i*lnxstep)),log(1-exp(log(1e-6)+(i+1)*lnxstep))};
			g0[0]=log(gridA[(xMax*j+i)*7+Iparton]);
			g0[1]=log(gridA[(xMax*j+i+1)*7+Iparton]);
			j++;
			g1[0]=log(gridA[(xMax*j+i)*7+Iparton]);
			g1[1]=log(gridA[(xMax*j+i+1)*7+Iparton]);
			j++;
			g2[0]=log(gridA[(xMax*j+i)*7+Iparton]);
			g2[1]=log(gridA[(xMax*j+i+1)*7+Iparton]);
			g[0]=exp(fitLinear(log(1-x),vln1_x,g0));
			g[1]=exp(fitLinear(log(1-x),vln1_x,g1));
			g[2]=exp(fitLinear(log(1-x),vln1_x,g2));
		}
		//if x<1e-5, we use A*x^B to do the interpolation
		//for valance quark, B>0; for gluon and sea quark, B<0
		else if(x<1e-4)
		{
			double vlnx[2]={(double)i,(double)(i+1)};
			g0[0]=log(gridA[(xMax*j+i)*7+Iparton]);
			g0[1]=log(gridA[(xMax*j+i+1)*7+Iparton]);
			j++;
			g1[0]=log(gridA[(xMax*j+i)*7+Iparton]);
			g1[1]=log(gridA[(xMax*j+i+1)*7+Iparton]);
			j++;
			g2[0]=log(gridA[(xMax*j+i)*7+Iparton]);
			g2[1]=log(gridA[(xMax*j+i+1)*7+Iparton]);
			g[0]=exp(fitLinear(lnx,vlnx,g0));
			g[1]=exp(fitLinear(lnx,vlnx,g1));
			g[2]=exp(fitLinear(lnx,vlnx,g2));
		}
		//we use quadratic interpolation method for other situations
		else
		{
			double vlnx[3]={(double)i,(double)(i+1),(double)(i+2)};
			g0[0]=gridA[(xMax*j+i)*7+Iparton];
			g0[1]=gridA[(xMax*j+i+1)*7+Iparton];
			g0[2]=gridA[(xMax*j+i+2)*7+Iparton];
			j++;
			g1[0]=gridA[(xMax*j+i)*7+Iparton];
			g1[1]=gridA[(xMax*j+i+1)*7+Iparton];
			g1[2]=gridA[(xMax*j+i+2)*7+Iparton];
			j++;
			g2[0]=gridA[(xMax*j+i)*7+Iparton];
			g2[1]=gridA[(xMax*j+i+1)*7+Iparton];
			g2[2]=gridA[(xMax*j+i+2)*7+Iparton];
			g[0]=fitQuadratic(lnx,vlnx,g0);
			g[1]=fitQuadratic(lnx,vlnx,g1);
			g[2]=fitQuadratic(lnx,vlnx,g2);
		}
		//if Q2>1, we do the interpolation to the variable ln(Q^2)
		if(Q2>1)
		{
			double vlnQ2[3]={(double)(j-2),(double)(j-1),(double)j};
			return fitQuadratic(lnQ2,vlnQ2,g);
		}
		//if Q2<1, we do the interpolation to the variable Q^2
		else 
		{
			double vQ2[3]={0.125*pow(2,j-2),0.125*pow(2,j-1),0.125*pow(2,j)};
			return fitQuadratic(Q2,vQ2,g);
		}	
	}
}

//return the parton distributions of deuteron of different kinds at x and Q^2
double nIMParton::getPDFD(int Iparton, double x, double Q2) const
{
        if(Iparton==-3 || Iparton==3)return getXSSeaD(x,Q2)/2.0/x;
        else if(Iparton==-2)return getXDSeaD(x,Q2)/2.0/x;
        else if(Iparton==2)return getXDSeaD(x,Q2)/2.0/x+getXDVD(x,Q2)/x;
        else if(Iparton==-1)return getXUSeaD(x,Q2)/2.0/x;
        else if(Iparton==1)return getXUSeaD(x,Q2)/2.0/x+getXUVD(x,Q2)/x;
        else if(Iparton==0)return getXGluonD(x,Q2)/x;
        else
        {
                cout<<"!!->Error: Unknown Iparton type."<<" (Iparton = "<<Iparton<<"?)"<<endl;
                cout<<"!!->Iparton should be one of these: [-3,-2, ... 2, 3]"<<endl;
                return 0;
        }

}

//a method which returns xuv of deuteron
double nIMParton::getXUVD(double x, double Q2) const
{
	return getPDFTypeD(1,x,Q2); 
}

//a method which returns xdv of deuteron
double nIMParton::getXDVD(double x, double Q2) const
{
        return getPDFTypeD(2,x,Q2);
}

//a method which returns xusea of deuteron
double nIMParton::getXUSeaD(double x, double Q2) const
{
        return getPDFTypeD(3,x,Q2);
}

//a method which returns xdsea of deuteron
double nIMParton::getXDSeaD(double x, double Q2) const
{
        return getPDFTypeD(4,x,Q2);
}

//a method which returns xssea of deuteron
double nIMParton::getXSSeaD(double x, double Q2) const
{
        return getPDFTypeD(5,x,Q2);
}

//a method which returns structure function F2 of deuteron
double nIMParton::getF2D(double x, double Q2) const
{
        return getPDFTypeD(6,x,Q2);
}

//a method which returns xgluon of deuteron
double nIMParton::getXGluonD(double x, double Q2) const
{
        return getPDFTypeD(0,x,Q2);
}

//a method which returns different types of distributions of deuteron
double nIMParton::getPDFTypeD(int Iparton, double x, double Q2) const
{
	if(Iparton<0 || Iparton>6)
	{
		cout<<"!!->Error: Wrong Iparton input for getPDFType(int Iparton, double x, double Q2)."<<endl;
		return 0;
	}
	else
	{
		double lnx, lnQ2;
		int i=(int)(lnx=log(x*1e6)/lnxstep);
		int j=(int)(lnQ2=log(Q2*8)/lnQ2step);
		double g0[3], g1[3], g2[3], g[3]={0};
		if(i<0)i=0;
		if(i>(int)(xMax-3))i=xMax-3;
		if(j<0)j=0;
		if(j>29)j=29;
		//avoid log(1-x) calculation in below algorithm
		if(x>0.9999)return 0.0;
		//if x>0.5, we use A(1-x)^B to do the interpolation
		else if(x>0.5)
		{
			double vln1_x[2]={log(1-exp(log(1e-6)+i*lnxstep)),log(1-exp(log(1e-6)+(i+1)*lnxstep))};
			g0[0]=log(gridD[(xMax*j+i)*7+Iparton]);
			g0[1]=log(gridD[(xMax*j+i+1)*7+Iparton]);
			j++;
			g1[0]=log(gridD[(xMax*j+i)*7+Iparton]);
			g1[1]=log(gridD[(xMax*j+i+1)*7+Iparton]);
			j++;
			g2[0]=log(gridD[(xMax*j+i)*7+Iparton]);
			g2[1]=log(gridD[(xMax*j+i+1)*7+Iparton]);
			g[0]=exp(fitLinear(log(1-x),vln1_x,g0));
			g[1]=exp(fitLinear(log(1-x),vln1_x,g1));
			g[2]=exp(fitLinear(log(1-x),vln1_x,g2));
		}
		//if x<1e-5, we use A*x^B to do the interpolation
		//for valance quark, B>0; for gluon and sea quark, B<0
		else if(x<1e-4)
		{
			double vlnx[2]={(double)i,(double)(i+1)};
			g0[0]=log(gridD[(xMax*j+i)*7+Iparton]);
			g0[1]=log(gridD[(xMax*j+i+1)*7+Iparton]);
			j++;
			g1[0]=log(gridD[(xMax*j+i)*7+Iparton]);
			g1[1]=log(gridD[(xMax*j+i+1)*7+Iparton]);
			j++;
			g2[0]=log(gridD[(xMax*j+i)*7+Iparton]);
			g2[1]=log(gridD[(xMax*j+i+1)*7+Iparton]);
			g[0]=exp(fitLinear(lnx,vlnx,g0));
			g[1]=exp(fitLinear(lnx,vlnx,g1));
			g[2]=exp(fitLinear(lnx,vlnx,g2));
		}
		//we use quadratic interpolation method for other situations
		else
		{
			double vlnx[3]={(double)i,(double)(i+1),(double)(i+2)};
			g0[0]=gridD[(xMax*j+i)*7+Iparton];
			g0[1]=gridD[(xMax*j+i+1)*7+Iparton];
			g0[2]=gridD[(xMax*j+i+2)*7+Iparton];
			j++;
			g1[0]=gridD[(xMax*j+i)*7+Iparton];
			g1[1]=gridD[(xMax*j+i+1)*7+Iparton];
			g1[2]=gridD[(xMax*j+i+2)*7+Iparton];
			j++;
			g2[0]=gridD[(xMax*j+i)*7+Iparton];
			g2[1]=gridD[(xMax*j+i+1)*7+Iparton];
			g2[2]=gridD[(xMax*j+i+2)*7+Iparton];
			g[0]=fitQuadratic(lnx,vlnx,g0);
			g[1]=fitQuadratic(lnx,vlnx,g1);
			g[2]=fitQuadratic(lnx,vlnx,g2);
		}
		//if Q2>1, we do the interpolation to the variable ln(Q^2)
		if(Q2>1)
		{
			double vlnQ2[3]={(double)(j-2),(double)(j-1),(double)j};
			return fitQuadratic(lnQ2,vlnQ2,g);
		}
		//if Q2<1, we do the interpolation to the variable Q^2
		else 
		{
			double vQ2[3]={0.125*pow(2,j-2),0.125*pow(2,j-1),0.125*pow(2,j)};
			return fitQuadratic(Q2,vQ2,g);
		}	
	}
}

//return the parton distributions of free nucleon of different kinds at x and Q^2
double nIMParton::getPDFN(int Iparton, double x, double Q2) const
{
        if(Iparton==-3 || Iparton==3)return getXSSeaN(x,Q2)/2.0/x;
        else if(Iparton==-2)return getXDSeaN(x,Q2)/2.0/x;
        else if(Iparton==2)return getXDSeaN(x,Q2)/2.0/x+getXDVN(x,Q2)/x;
        else if(Iparton==-1)return getXUSeaN(x,Q2)/2.0/x;
        else if(Iparton==1)return getXUSeaN(x,Q2)/2.0/x+getXUVN(x,Q2)/x;
        else if(Iparton==0)return getXGluonN(x,Q2)/x;
        else
        {
                cout<<"!!->Error: Unknown Iparton type."<<" (Iparton = "<<Iparton<<"?)"<<endl;
                cout<<"!!->Iparton should be one of these: [-3,-2, ... 2, 3]"<<endl;
                return 0;
        }

}

//a method which returns xuv of free nucleon
double nIMParton::getXUVN(double x, double Q2) const
{
	return getPDFTypeN(1,x,Q2); 
}

//a method which returns xdv of free nucleon
double nIMParton::getXDVN(double x, double Q2) const
{
        return getPDFTypeN(2,x,Q2);
}

//a method which returns xusea of free nucleon
double nIMParton::getXUSeaN(double x, double Q2) const
{
        return getPDFTypeN(3,x,Q2);
}

//a method which returns xdsea of free nucleon
double nIMParton::getXDSeaN(double x, double Q2) const
{
        return getPDFTypeN(4,x,Q2);
}

//a method which returns xssea of free nucleon
double nIMParton::getXSSeaN(double x, double Q2) const
{
        return getPDFTypeN(5,x,Q2);
}

//a method which returns structure function F2 of free nucleon
double nIMParton::getF2N(double x, double Q2) const
{
        return getPDFTypeN(6,x,Q2);
}

//a method which returns xgluon of free nucleon
double nIMParton::getXGluonN(double x, double Q2) const
{
        return getPDFTypeN(0,x,Q2);
}

//a method which returns different types of distributions of free nucleon
double nIMParton::getPDFTypeN(int Iparton, double x, double Q2) const
{
	if(Iparton<0 || Iparton>6)
	{
		cout<<"!!->Error: Wrong Iparton input for getPDFType(int Iparton, double x, double Q2)."<<endl;
		return 0;
	}
	else
	{
		double lnx, lnQ2;
		int i=(int)(lnx=log(x*1e6)/lnxstep);
		int j=(int)(lnQ2=log(Q2*8)/lnQ2step);
		double g0[3], g1[3], g2[3], g[3]={0};
		if(i<0)i=0;
		if(i>(int)(xMax-3))i=xMax-3;
		if(j<0)j=0;
		if(j>29)j=29;
		//avoid log(1-x) calculation in below algorithm
		if(x>0.9999)return 0.0;
		//if x>0.5, we use A(1-x)^B to do the interpolation
		else if(x>0.5)
		{
			double vln1_x[2]={log(1-exp(log(1e-6)+i*lnxstep)),log(1-exp(log(1e-6)+(i+1)*lnxstep))};
			g0[0]=log(gridN[(xMax*j+i)*7+Iparton]);
			g0[1]=log(gridN[(xMax*j+i+1)*7+Iparton]);
			j++;
			g1[0]=log(gridN[(xMax*j+i)*7+Iparton]);
			g1[1]=log(gridN[(xMax*j+i+1)*7+Iparton]);
			j++;
			g2[0]=log(gridN[(xMax*j+i)*7+Iparton]);
			g2[1]=log(gridN[(xMax*j+i+1)*7+Iparton]);
			g[0]=exp(fitLinear(log(1-x),vln1_x,g0));
			g[1]=exp(fitLinear(log(1-x),vln1_x,g1));
			g[2]=exp(fitLinear(log(1-x),vln1_x,g2));
		}
		//if x<1e-5, we use A*x^B to do the interpolation
		//for valance quark, B>0; for gluon and sea quark, B<0
		else if(x<1e-4)
		{
			double vlnx[2]={(double)i,(double)(i+1)};
			g0[0]=log(gridN[(xMax*j+i)*7+Iparton]);
			g0[1]=log(gridN[(xMax*j+i+1)*7+Iparton]);
			j++;
			g1[0]=log(gridN[(xMax*j+i)*7+Iparton]);
			g1[1]=log(gridN[(xMax*j+i+1)*7+Iparton]);
			j++;
			g2[0]=log(gridN[(xMax*j+i)*7+Iparton]);
			g2[1]=log(gridN[(xMax*j+i+1)*7+Iparton]);
			g[0]=exp(fitLinear(lnx,vlnx,g0));
			g[1]=exp(fitLinear(lnx,vlnx,g1));
			g[2]=exp(fitLinear(lnx,vlnx,g2));
		}
		//we use quadratic interpolation method for other situations
		else
		{
			double vlnx[3]={(double)i,(double)(i+1),(double)(i+2)};
			g0[0]=gridN[(xMax*j+i)*7+Iparton];
			g0[1]=gridN[(xMax*j+i+1)*7+Iparton];
			g0[2]=gridN[(xMax*j+i+2)*7+Iparton];
			j++;
			g1[0]=gridN[(xMax*j+i)*7+Iparton];
			g1[1]=gridN[(xMax*j+i+1)*7+Iparton];
			g1[2]=gridN[(xMax*j+i+2)*7+Iparton];
			j++;
			g2[0]=gridN[(xMax*j+i)*7+Iparton];
			g2[1]=gridN[(xMax*j+i+1)*7+Iparton];
			g2[2]=gridN[(xMax*j+i+2)*7+Iparton];
			g[0]=fitQuadratic(lnx,vlnx,g0);
			g[1]=fitQuadratic(lnx,vlnx,g1);
			g[2]=fitQuadratic(lnx,vlnx,g2);
		}
		//if Q2>1, we do the interpolation to the variable ln(Q^2)
		if(Q2>1)
		{
			double vlnQ2[3]={(double)(j-2),(double)(j-1),(double)j};
			return fitQuadratic(lnQ2,vlnQ2,g);
		}
		//if Q2<1, we do the interpolation to the variable Q^2
		else 
		{
			double vQ2[3]={0.125*pow(2,j-2),0.125*pow(2,j-1),0.125*pow(2,j)};
			return fitQuadratic(Q2,vQ2,g);
		}	
	}
}


//linear interpolation method
double nIMParton::fitQuadratic(double x, double* px, double* pf) const
{
	double f01=(pf[1]-pf[0])/(px[1]-px[0]);
	double f12=(pf[2]-pf[1])/(px[2]-px[1]);
	double f012=(f12-f01)/(px[2]-px[0]);
	return pf[0]+f01*(x-px[0])+f012*(x-px[0])*(x-px[1]);
}

//quadratic interpolation method
double nIMParton::fitLinear(double x, double* px, double* pf) const
{
	double f01=(pf[1]-pf[0])/(px[1]-px[0]);
	return pf[0]+f01*(x-px[0]);
}



