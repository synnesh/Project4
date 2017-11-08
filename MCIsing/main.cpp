#include <iostream>

using namespace std;

void metropolis(int, int, int, int);
void initialize(int);
void output(int, int, double, double, double, double, double);
double ran2(long *);

int main(int argc, char* argv[])
{
    int latticeDim;int mc_cycles;
    double initial_temp;double max_temp;double temp_step;
    long idum;
    double w[17];


    idum=-1;

    double initial_temp = 2.0; double max_temp=2.4;double temp_step=0.1;
    mc_cycles=100000;
    latticeDim=2;
    for (int de=-8;de<=8;de+=4){
        w[de+8]=exp(-de/T);
    }


    //Loop over temperatures
    for(double T =initial_temp;T<=max_temp;T+=temp_step)
        //Precalculate possible changes in probability
        for (int de=-8;de<=8;de+=4){
            w[de+8]=exp(-de/T);
        }
        //Reset Energy and magnetization
        totalE=0;totalE2=0;
        totalM=0;totalM2=0;

        //Initialize random lattice
        intialize(latticeDim);
        //MC loop
        for (int i = 0; i<mc_cycles;i++){
            metropolis(latticeDim, E, M);
            totalE+=E; totalE2+=E*E;
            totalM+=M; totalM2+=M*M; m_abs +=fabs(M);
        }
        output(latticeDim, i, T, totalE, totalM, m_abs, totalE2, totalM2);



}


void metropolis(int latticeDim, int E, int M, int idum, int spinmatrix){
    for (int j =0;j<latticeDim;j++){
        for(int k=0;k<latticeDim;j++){
            x= (int) ran2(*idum)*(double)latticeDim;
            y= (int) ran2(*idum)*(double)latticeDim;

            deltaE=2*spinmatrix[x][y](spinmatrix[x][periodic(y, latticeDim,1)]+spinmatrix[periodic(x, latticeDim,1)][y]+spinmatrix[x][periodic(y, latticeDim, -1)]+spinmatrix[periodic(x, latticeDim,-1][y]);

            if(ran2(idum) <= w[deltaE+8]){
                spinmatrix[x][y]*=-1;
                M+=2*spinmatrix[x][y];
                E-=(double) deltaE;
            }
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
