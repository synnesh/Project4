#include <iostream>
#include "omp.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <armadillo>

using namespace arma;
using namespace std;

ofstream ofile;

inline int periodic(int location, int latticeDim,int move){
    return (location+latticeDim+move)%latticeDim;
}


void metropolis(int, double&, double&, long *, mat, vec);
void initialize(int, long *, double&, double&, mat&);
void output(int, int, double, double, double, double, double, double);
double ran2(long *);

int main(int argc, char* argv[])
{
    char* outfilename;
    int latticeDim;int mc_cycles;
    double initial_temp;double max_temp;double temp_step;
    long idum;

    outfilename= argv[1];
    ofile.open(outfilename);

    idum=-1;

    initial_temp = 2.0; max_temp=2.1;temp_step=0.1;
    mc_cycles=10;
    latticeDim=2;

    double totalE,totalE2,totalM, totalM2, m_abs, E,M;

    //Loop over temperatures
    vec w(17);
    for (int i =0;i<17;i++){
        w[i]=0;
    }

    mat spinmatrix(latticeDim, latticeDim);
    for(double T =initial_temp;T<=max_temp;T+=temp_step){
        //Precalculate possible changes in probability
        for (int de=-8;de<=8;de+=4){
            w[de+8]=exp(-de/T);
        }
        for (int i=0;i<17;i++){
            cout<<w[i] <<" ";
        }
        cout <<"\n";
        //Reset Energy and magnetization
        totalE=0;totalE2=0;
        totalM=0;totalM2=0;
        E=M=0;
        //Initialize random lattice
        initialize(latticeDim, &idum, E, M, spinmatrix);
        //MC loop

        for (int i = 1; i<mc_cycles;i++){

            metropolis(latticeDim, E, M, &idum, spinmatrix, w);
            totalE+=E; totalE2+=E*E;
            totalM+=M; totalM2+=M*M; m_abs +=fabs(M);
            //Write variables to file
            output(latticeDim, i, T, totalE, totalM, m_abs, totalE2, totalM2);
        }

    }

return 0;
}


void metropolis(int latticeDim, double& E, double& M, long *idum, mat spinmatrix, vec w){
    for (int j =0;j<latticeDim;j++){
        for(int k=0;k<latticeDim;k++){

            int x= (int) ran2(idum)*(double)latticeDim;
            int y= (int) ran2(idum)*(double)latticeDim;


            double deltaE=2*spinmatrix(x,y)*(spinmatrix(x,periodic(y, latticeDim,1))+spinmatrix(periodic(x, latticeDim,1),y)+spinmatrix(x,periodic(y, latticeDim, -1))+spinmatrix(periodic(x, latticeDim,-1),y));
            cout<<deltaE<<"\n";
            if(ran2(idum) <= w[deltaE+8]){
                spinmatrix(x,y)*=-1;
                M+=2*spinmatrix(x,y);
                E+= deltaE;
            }
            for(int i=0;i<latticeDim; i++){
                for(int k=0;k<latticeDim; k++){
                    cout <<spinmatrix(i,k)<<"  ";
                }
                cout<<"\n";
            }
            cout<<"\n";
        }

    }

}

void initialize(int latticeDim, long *idum, double& E, double &M, mat& spinmatrix)
{
    //Construct matrix with random spin configuration

    for (int x=0;x<latticeDim;x++){
        for(int y=0;y<latticeDim;y++){
            double s=ran2(idum);
            if(s<.5){
                spinmatrix(x,y)=1.;
            }
            else{
                spinmatrix(x,y)=-1.;
            }
        }
    }
    //Calculate E and M

    for (int x=0;x<latticeDim;x++){
        for (int y = 0;y<latticeDim;y++){
               M+=spinmatrix(x,y);
               E+=spinmatrix(x,y)*spinmatrix(x, periodic(y,latticeDim,1))+spinmatrix(x,y)*spinmatrix(periodic(x, latticeDim, 1));
        }

    }
}

void output(int n_spins, int mcs, double temperature, double totalE, double totalM, double Mabs, double totalE2, double totalM2)
{
  double norm = 1/((double) (mcs));  // divided by total number of cycles
  double Etotal_average = totalE*norm;
  double E2total_average = totalE2*norm;
  double Mtotal_average = totalM*norm;
  double M2total_average = totalM2*norm;
  double Mabstotal_average = Mabs*norm;
  // all expectation values are per spin, divide by 1/n_spins/n_spins
  double Evariance = (E2total_average- Etotal_average*Etotal_average)/n_spins/n_spins;
  double Mvariance = (M2total_average - Mabstotal_average*Mabstotal_average)/n_spins/n_spins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << temperature;
  ofile << setw(15) << setprecision(8) << Etotal_average/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
  ofile << setw(15) << setprecision(8) << Mtotal_average/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Mvariance/temperature;
  ofile << setw(15) << setprecision(8) << Mabstotal_average/n_spins/n_spins << endl;
} // end output function



#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
  int            j;
  long           k;
  static long    idum2 = 123456789;
  static long    iy=0;
  static long    iv[NTAB];
  double         temp;

  if(*idum <= 0) {
    if(-(*idum) < 1) *idum = 1;
    else             *idum = -(*idum);
    idum2 = (*idum);
    for(j = NTAB + 7; j >= 0; j--) {
      k     = (*idum)/IQ1;
      *idum = IA1*(*idum - k*IQ1) - k*IR1;
      if(*idum < 0) *idum +=  IM1;
      if(j < NTAB)  iv[j]  = *idum;
    }
    iy=iv[0];
  }
  k     = (*idum)/IQ1;
  *idum = IA1*(*idum - k*IQ1) - k*IR1;
  if(*idum < 0) *idum += IM1;
  k     = idum2/IQ2;
  idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
  if(idum2 < 0) idum2 += IM2;
  j     = iy/NDIV;
  iy    = iv[j] - idum2;
  iv[j] = *idum;
  if(iy < 1) iy += IMM1;
  if((temp = AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
