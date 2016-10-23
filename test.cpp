#include<iostream>
#include<fstream>
extern "C"{
#include<math.h>
#include<time.h>
#include<unistd.h>
#include<sys/time.h>
}
#include "nIMParton.h"
using namespace std;


//nIMParton He3(2,3);                // Z=2  and A=3,   which is ^3Helium
//nIMParton Al27(13,27);             // Z=13 and A=27,  which is ^27Aluminum
//nIMParton Ca40(20,40);             // Z=20 and A=40,  which is ^40Calcium
//nIMParton Fe56(26,56);             // Z=26 and A=56,  which is ^56Iron
nIMParton Pb208(82,208);           // Z=82 and A=208, which is ^208Lead

int main()
{
    struct timeval tstart,tend;
    gettimeofday(&tstart,NULL);                 //get the starting time

    double lnxstep=(log(1)-log(0.1))/100;       //step length of log(x) for output data
    ofstream outdata1("GluonRatio_Q2_10.dat");  //output file for nuclear modification of gluon distribution
    ofstream outdata2("uRatio_Q2_10.dat");      //output file for nuclear modification of up quark distribution
    ofstream outdata3("ubarRatio_Q2_10.dat");   //output file for nuclear modification of ubar distribution
    ofstream outdata4("dRatio_Q2_10.dat");      //output file for nuclear modification of down quark distribution
    ofstream outdata5("dbarRatio_Q2_10.dat");   //output file for nuclear modification of dbar distribution
    ofstream outdata6("sRatio_Q2_10.dat");      //output file for nuclear modification of s distribuiton
    ofstream outdata7("F2Ratio_Q2_10.dat");     //output file for structure function F2 ratios

    //using data set A
    Pb208.setDataSet(1);             //1 for data set A and 2 data set B
    //get the nuclear modifications of lead target
    for(int i=0;i<600;i++)
    {
        double x = exp(log(1e-6)+i*lnxstep);
        outdata1<<x<<" "<<Pb208.getRToD(0,  x, 10)<<endl;
        outdata2<<x<<" "<<Pb208.getRToD(1,  x, 10)<<endl;
        outdata3<<x<<" "<<Pb208.getRToD(-1, x, 10)<<endl;
        outdata4<<x<<" "<<Pb208.getRToD(2,  x, 10)<<endl;
        outdata5<<x<<" "<<Pb208.getRToD(-2, x, 10)<<endl;
        outdata6<<x<<" "<<Pb208.getRToD(3,  x, 10)<<endl;
        outdata7<<x<<" "<<Pb208.getF2RToD(x, 10)<<endl;
    }

    gettimeofday(&tend,NULL);         //get the ending time
    double timeuse=1000000*(tend.tv_sec-tstart.tv_sec)+(tend.tv_usec-tstart.tv_usec);
    //display the total program running time
    cout<<"    Runtime : "<<(timeuse/1000.0)<<" ms."<<endl;
}
