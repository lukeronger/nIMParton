#include<iostream>
#include<fstream>
extern "C"{
#include<math.h>
#include<time.h>
#include<unistd.h>
#include<sys/time.h>
}

//// apply nIMParton class
#include "nIMParton.h"

using namespace std;


int main()
{
	// construct a nIMParton object
	nIMParton He3(2,3);                         // Z=2 and A=3, which is ^3Helium

	struct timeval tstart, tend;
	gettimeofday(&tstart,NULL);                 //get the starting time

	double lnxstep=(log(1)-log(0.1))/100;       //step length of log(x) for output data table
	ofstream outdata1("GluonRatio_Q2_10_ToN_SetB.dat");  //output file for gluon ratio of nucleus to free nucleon
	ofstream outdata2("uRatio_Q2_10_ToN_SetB.dat");      //output file for up quark ratio of nucleus to free nucleon
	ofstream outdata3("ubarRatio_Q2_10_ToN_SetB.dat");   //output file for ubar ratio of nucleus to free nucleon
	ofstream outdata4("GluonRatio_Q2_10_ToD_SetB.dat");  //output file for gluon ratio of nucleus to Deuterium
	ofstream outdata5("uRatio_Q2_10_ToD_SetB.dat");      //output file for up quark ratio of nucleus to Deuterium
	ofstream outdata6("ubarRatio_Q2_10_ToD_SetB.dat");   //output file for ubar ratio of nucleus to Deuterium

	//using data set B
	He3.setDataSet(2);                          //1 for data set A and 2 for data set B
	//get the nuclear modifications of ^3He target
        //at Q^2 = 10 GeV^2
	for(int i=0;i<600;i++)
	{
		double x = exp(log(1e-6)+i*lnxstep);
		//get gluon ratio of the bound proton in ^3He to free proton
		outdata1<<x<<"\t"<<He3.getRToN_p(0,  x, 10)<<endl;
		//get up quark ratio of the bound proton in ^3He to free proton
		outdata2<<x<<"\t"<<He3.getRToN_p(1,  x, 10)<<endl;
		//get ubar ratio of the bound proton in ^3He to free proton
		outdata3<<x<<"\t"<<He3.getRToN_p(-1, x, 10)<<endl;
		//get gluon ratio of the bound proton in ^3He to the bound proton in Deuterium
		outdata4<<x<<"\t"<<He3.getRToD_p(0,  x, 10)<<endl;
		//get up quark ratio of the bound proton in ^3He to the bound proton in Deuterium
		outdata5<<x<<"\t"<<He3.getRToD_p(1,  x, 10)<<endl;
		//get ubar ratio of the bound proton in ^3He to the bound proton in Deuterium
		outdata6<<x<<"\t"<<He3.getRToD_p(-1, x, 10)<<endl;
	}

	gettimeofday(&tend,NULL);                   //get the ending time
	double timeuse = 1000000*(tend.tv_sec-tstart.tv_sec) + (tend.tv_usec-tstart.tv_usec);
	//display the total time for calling the nIMParton member function
	cout<<"    Runtime : "<<(timeuse/1000.0)<<" ms."<<endl;
	
	return 0;
}
