#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <ctime>
#include <chrono>
#include <time.h>

#if(INTEL==0)
	#include <boost/filesystem.hpp>
#endif

#define PI (3.141592653589793)
#define RKON 1

//*********************************************************************************
//	Includes for FFTW
//*********************************************************************************

#include <fftw3.h>

//*********************************************************************************
//	Includes for the random number generators
//*********************************************************************************

#include "3_RNG_IMPLEMENTATION.c"


#define TESTING (0)

int main(int argc,char *argv[]){


    int64_t TAULEN=atoi(argv[1]);		   //NUMBER OF STEPS ALONG THE EUCLIDEAN TIME DIRECTION
    int64_t POSLEN=atoi(argv[18]);        //steps in eucl x-direction

    int64_t OMGLEN=atoi(argv[2]);		   //NUMBER OF STEPS ALONG THE IMAGIANRY FREQUENCY DIRECTION
    int64_t MOMLEN=atoi(argv[19]);        //steps in imag x-momentum

    int64_t INTERACTION=atoi(argv[3]);          //WHETHER PHI^4 TERM IS INCLUDED OR NOT

    double BETA=atof(argv[4]);			//INVERSE TEMPERATURE OF THE SYSTEM
    double LENX=atof(argv[20]);         //dimensions of the box

    double DT=BETA/(double)TAULEN;                    //LATTICE SPACING IN EUCLIDEAN TIME
    double DX=LENX/(double)POSLEN;                    //deltax

    double MASS=atof(argv[5]);			//MASS OF THE REAL SCALAR FIELD
    double LAMBDA=atof(argv[6]);		//SELF COUPLING OF THE REAL SCALAR FIELD

    int64_t BINLEN=atof(argv[21]);              //LENGHT OF BINNING VECTOR

    int64_t VERBOSE=atoi(argv[14]);             //WHETHER TO OUTPUT DEBUGGING INFORMATION

    double dt5=atof(argv[7]);            //STEP SIZE FOR THE STOCHASTIC QUANTIZATION LANGEVIN

    double STEPS=atof(argv[8]);                 //NUMBER OF STEPS TO PROCEED IN THE LANGEVIN EVOLUTION

	int64_t METSTEPS=STEPS;
    int64_t STARTMEAS=1/dt5;			//AFTER HOW MANY STEPS TO START THE MEASUREMENTS
    int64_t MEASEACH=atoi(argv[9]);		//MEASURE THE OBSERVABLES AT EACH MEASEACH STEPS
    int64_t MEASCORREACH=atoi(argv[10]);	//MEASURE THE OBSERVABLES AT EACH MEASEACH STEPS
    int64_t MEASBAYESCORREACH=atoi(argv[11]);

	int64_t PRINTCURRENT=atof(argv[23]);

	//----------------------
	//FIELDS AND CORRELATORS
	//----------------------
	double* FIELDcurr = new double [TAULEN*POSLEN*POSLEN];
	double* FIELDnext = new double [TAULEN*POSLEN*POSLEN];
	double* FIELDhalf = new double [TAULEN*POSLEN*POSLEN];
	double* NOISE = new double [TAULEN*POSLEN*POSLEN];
	double* Drift = new double [TAULEN*POSLEN*POSLEN];

	double* EuclCORR = new double [TAULEN*POSLEN*POSLEN];
	double* EuclCORRAvg = new double [TAULEN*POSLEN*POSLEN];
	double* EuclCORRZero = new double [TAULEN];
	double* EuclCORRZeroAvg = new double [TAULEN];

	double* FIELDOMGA = new double [2*OMGLEN*MOMLEN*MOMLEN];
	double* FIELDOMGAINT = new double [2*OMGLEN*MOMLEN*MOMLEN];
	double* FIELDOMGAhalf = new double [2*OMGLEN*MOMLEN*MOMLEN];
	double* FOURNOISE = new double [2*OMGLEN*MOMLEN*MOMLEN];
	double* IMFreqDrift = new double [2*OMGLEN*MOMLEN*MOMLEN];

	double* ReCORR = new double [OMGLEN*MOMLEN*MOMLEN];
	double* ImCORR = new double [OMGLEN*MOMLEN*MOMLEN];
	double* ReCORRmom = new double [OMGLEN*BINLEN];
	double* ImCORRmom = new double [OMGLEN*BINLEN];
	double* ReCORRmomAvg = new double [OMGLEN*BINLEN];
	double* ImCORRmomAvg = new double [OMGLEN*BINLEN];

	//---------------
	//MOMENTUM ARRAYS
	//---------------
	double* momtem = new double [OMGLEN];
	double* momx = new double [MOMLEN];
	double* momsquared = new double [OMGLEN*MOMLEN*MOMLEN];
	double* momxsquared = new double [MOMLEN*MOMLEN];
	double* momval = new double [BINLEN];
	int64_t* mombin = new int64_t [BINLEN];


	//-----------
	//MISCELLANEA
	//-----------
	double mommax = 0.;
	double DBIN;
	int64_t PCNT,
			NCNT,
			negx,
			negy,
			bindex,
			id;
	double _twopi_over_beta=2.*PI/TAULEN;
	double conv=(double)TAULEN/(double)OMGLEN;
	double quadtrm = LAMBDA/24.;
	double h=0;
	char* fname = new char [256];

	//----------------------
	//SETUP RANDOM VARIABLES
	//----------------------
	double SEED=25775447298;
	TLP_ALGORITHM *GSLRAN;
	GSLRAN = (TLP_ALGORITHM *) malloc(sizeof(TLP_ALGORITHM));
	TLP_ALGORITHM_SET (GSLRAN, 2);
	TLP_RANDGEN *RNG;
	RNG = TLP_RANDGEN_INITIALIZE( GSLRAN, SEED );


	//------------------
	//INITIAL CONDITIONS
	//------------------
	for(int64_t i=0;i<BINLEN;i++){
		momval[i]=0.;
		mombin[i]=0;
	}
	for(int64_t i=0;i<OMGLEN*MOMLEN*MOMLEN;i++){
		ReCORR[i]=0.;
		ImCORR[i]=0.;
	}
	for(int64_t i=0;i<2*OMGLEN*MOMLEN*MOMLEN;i++){
		FIELDOMGA[i]=0.;
		FIELDOMGAhalf[i]=0.;
		FIELDOMGAINT[i]=0.;
		FOURNOISE[i]=0.;
	}
	for(int64_t i=0;i<TAULEN*POSLEN*POSLEN;i++){
		FIELDcurr[i]=0.;
		FIELDnext[i]=0.;
		FIELDhalf[i]=0.;
		EuclCORR[i]=0.;
		EuclCORRAvg[i]=0.;
		NOISE[i]=0.;
	}
	for(int64_t i=0;i<TAULEN;i++){
		EuclCORRZero[i]=0;
		EuclCORRZeroAvg[i]=0;
	}
	for(int64_t i=0;i<BINLEN*OMGLEN;i++){
	  ReCORRmomAvg[i]=0.;
	  ImCORRmomAvg[i]=0.;
	}

	//CHECK IF INTERACTION==1
	if(INTERACTION!=1 && INTERACTION!=0){
		std::cerr << "Set INTERACTION variable either 0 or 1. Abort." << std::endl;
		exit(-1);
	}

	//-----------------------
	//FOURIER TRANSFORM PLANS
	//-----------------------
	fftw_plan FTForInit,IFTForInteraction,FTForInteraction,FTForNoise,HalfFTForInteraction,FTEucl;

	FTForInit = fftw_plan_dft_3d(MOMLEN,MOMLEN,OMGLEN, (fftw_complex*)FIELDOMGA, (fftw_complex*)FIELDOMGA, FFTW_FORWARD, FFTW_ESTIMATE);
	IFTForInteraction = fftw_plan_dft_3d(MOMLEN,MOMLEN,OMGLEN, (fftw_complex*)FIELDOMGA, (fftw_complex*)FIELDOMGAINT, FFTW_BACKWARD, FFTW_ESTIMATE);
	HalfFTForInteraction = fftw_plan_dft_3d(MOMLEN,MOMLEN,OMGLEN, (fftw_complex*)FIELDOMGAhalf, (fftw_complex*)FIELDOMGAINT, FFTW_FORWARD, FFTW_ESTIMATE);
	FTForInteraction = fftw_plan_dft_3d(MOMLEN,MOMLEN,OMGLEN, (fftw_complex*)FIELDOMGAINT, (fftw_complex*)FIELDOMGAINT, FFTW_FORWARD, FFTW_ESTIMATE);
	FTForNoise = fftw_plan_dft_3d(MOMLEN,MOMLEN,OMGLEN, (fftw_complex*)FOURNOISE, (fftw_complex*)FOURNOISE,FFTW_FORWARD,FFTW_ESTIMATE);

	//------------------
	//COMPUTING MOMENTUM
	//------------------
	for(int64_t k=0;k<OMGLEN;k++) momtem[k]=sqrt(2.-2.*cos(k*conv*_twopi_over_beta))/(DT); //this is just k*2pi/OMGLEN
	for(int64_t k=0;k<MOMLEN;k++) momx[k]=sqrt(2.-2.*cos(k*2*PI/MOMLEN))/DX;

	//FULL MOMENTUM SQUARED
	for(int64_t x=0;x<MOMLEN;x++)
	for(int64_t y=0;y<MOMLEN;y++)
		for(int64_t k=0;k<OMGLEN;k++)
		  momsquared[(x*MOMLEN + y)*OMGLEN + k]=momtem[k]*momtem[k] + momx[x]*momx[x] + momx[y]*momx[y];


	//SPATIAL MOMENTUM SQUARED AND MAX MOMENTUM
	for(int64_t x=0;x<MOMLEN;x++)
		for(int64_t y=0;y<MOMLEN;y++){

			momxsquared[x*MOMLEN + y]=momx[x]*momx[x]+momx[y]*momx[y];
			if(momxsquared[x*MOMLEN + y] > mommax)
				mommax=momxsquared[x*MOMLEN + y];
	 	}

	DBIN=sqrt(mommax)/(double)BINLEN;

	//MOMENTUM BINNING
	for(int64_t k=0;k<MOMLEN*MOMLEN;k++){
		bindex = (int)floor(sqrt(momxsquared[k])/DBIN);

	if(bindex<BINLEN){

	  mombin[bindex]++;
	  momval[bindex] += sqrt(momxsquared[k])/DBIN;
	}
	else
	   std::cerr << k << std::endl;
	}

	//AVERAGING MOMVAL
	for(int64_t k=0;k<MOMLEN;k++)
		if(mombin[k]!=0)
			momval[k] /= mombin[k];

	//----------------------------
	//OUTPUT FOR BAYESIAN ANALYSIS
	//----------------------------
	std::fstream CORROUTE,CORROUTI;

	sprintf(fname,"%sEUCLCORRS.dat","./risultati");
    CORROUTE.open(fname, std::ios::out|std::ios::binary);
    sprintf(fname,"%sIMAGCORRS.dat","./risultati");
    CORROUTI.open(fname, std::ios::out|std::ios::binary);

	CORROUTE.write ((char*) &TAULEN, sizeof(int64_t));
	CORROUTI.write ((char*) &OMGLEN, sizeof(int64_t));

	for(int64_t k=0;k<TAULEN;k++) { double tau=DT*k; CORROUTE.write ((char*) &tau, sizeof(double)); }
	for(int64_t k=0;k<OMGLEN;k++) { CORROUTI.write ((char*) &momtem[k], sizeof(double)); }

	CORROUTE.flush();
	CORROUTI.flush();

	//--------------------------
	//BEGINNING OF LANGEVIN LOOP
	//--------------------------

	for(int64_t t5=0;t5<METSTEPS;t5++){

		if(t5%1000==0)
			fprintf(stderr, "StochQuant Step [%8ld]\r",t5);

		//GENERATING NOISE IN x-SPACE - for euclidean update
		for(int64_t i=0;i<TAULEN*POSLEN*POSLEN;i++)
			NOISE[i] = sqrt(2./DT)*(1./DX)*TLP_RAND_GAUSSIAN_ZIGGURAT (RNG, 1.);

		//GENERATING NOISE IN x-SPACE - for im-freq update
        for(int64_t w=0;w<OMGLEN*MOMLEN*MOMLEN;w++){
            FOURNOISE[2*w]=sqrt(2./DT)*(1./DX)*TLP_RAND_GAUSSIAN_ZIGGURAT (RNG, 1.);
            FOURNOISE[2*w+1]=0.;
        }

        //F-TRASFORM NOISE
        fftw_execute(FTForNoise);

        //NORMALIZING NOISE IN k-SPACE
        for(int64_t w=0;w<2*OMGLEN*MOMLEN*MOMLEN;w++)
            FOURNOISE[w]/=sqrt(OMGLEN*MOMLEN*MOMLEN);

//EUCLIDEAN UPDATE
		for(int64_t x=0;x<POSLEN;x++)
		   for(int64_t y=0;y<POSLEN;y++)
			   for(int64_t t=0;t<TAULEN;t++){

				   int64_t id =(x*POSLEN+y)*TAULEN + t;
				   double PhiDoubPrime = (1./(DT*DT))*(-2.*FIELDcurr[(x*POSLEN+y)*TAULEN + t] +                      //finite difference in t
														   FIELDcurr[(x*POSLEN+y)*TAULEN + (t+1)%TAULEN] +
														   FIELDcurr[(x*POSLEN+y)*TAULEN + (t-1+TAULEN)%TAULEN]) +

										 (1./(DX*DX))*(-2.*FIELDcurr[(x*POSLEN+y)*TAULEN + t] +                      //finite difference in x
														   FIELDcurr[(x*POSLEN+(y+1)%POSLEN)*TAULEN + t] +
														   FIELDcurr[(x*POSLEN+(y-1+POSLEN)%POSLEN)*TAULEN + t]) +

										 (1./(DX*DX))*(-2.*FIELDcurr[(x*POSLEN+y)*TAULEN + t] +                      //finite difference in y
														   FIELDcurr[(((x-1+POSLEN)%POSLEN)*POSLEN+y)*TAULEN + t] +
														   FIELDcurr[(((x+1)%POSLEN)*POSLEN+y)*TAULEN + t]);



				   Drift[id] = MASS*MASS*FIELDcurr[id] - PhiDoubPrime + INTERACTION*quadtrm*4*pow(FIELDcurr[id],3.);
				   FIELDnext[(x*POSLEN+y)*TAULEN + t] = FIELDcurr[(x*POSLEN+y)*TAULEN + t] - dt5*Drift[id] + sqrt(dt5)*NOISE[id];	//half Euler step
			   }

		//------------------------------------------------
        //COMPUTING EUCLIDEAN FIELD UPDATE - R-K 2nd order
        //------------------------------------------------
     	for(int64_t x=0;x<POSLEN;x++)
  			for(int64_t y=0;y<POSLEN;y++)
				for(int64_t t=0;t<TAULEN;t++){

					id = (x*POSLEN+y)*TAULEN + t;
					double PhiDoubPrime = (1./(DT*DT))*(-2.*FIELDcurr[(x*POSLEN+y)*TAULEN + t] +                      //finite difference in t
															FIELDcurr[(x*POSLEN+y)*TAULEN + (t+1)%TAULEN] +
															FIELDcurr[(x*POSLEN+y)*TAULEN + (t-1+TAULEN)%TAULEN]) +

										  (1./(DX*DX))*(-2.*FIELDcurr[(x*POSLEN+y)*TAULEN + t] +                      //finite difference in x
															FIELDcurr[(x*POSLEN+(y+1)%POSLEN)*TAULEN + t] +
															FIELDcurr[(x*POSLEN+(y-1+POSLEN)%POSLEN)*TAULEN + t]) +

										  (1./(DX*DX))*(-2.*FIELDcurr[(x*POSLEN+y)*TAULEN + t] +                      //finite difference in y
															FIELDcurr[(((x-1+POSLEN)%POSLEN)*POSLEN+y)*TAULEN + t] +
															FIELDcurr[(((x+1)%POSLEN)*POSLEN+y)*TAULEN + t]);

					Drift[id] = MASS*MASS*FIELDcurr[(x*POSLEN+y)*TAULEN + t] - PhiDoubPrime + INTERACTION*quadtrm*4*pow(FIELDcurr[(x*POSLEN+y)*TAULEN + t],3.);
					FIELDhalf[(x*POSLEN+y)*TAULEN + t] = FIELDcurr[(x*POSLEN+y)*TAULEN + t] - dt5*0.5*Drift[id] + sqrt(dt5*0.5)*NOISE[id];
				}

		for(int64_t x=0;x<POSLEN;x++)
  			for(int64_t y=0;y<POSLEN;y++)
				for(int64_t t=0;t<TAULEN;t++){

					id = (x*POSLEN+y)*TAULEN + t;
					double PhiDoubPrimeHalf = (1./(DT*DT))*(-2.*FIELDhalf[(x*POSLEN+y)*TAULEN + t] +                      //finite difference in t
																FIELDhalf[(x*POSLEN+y)*TAULEN + (t+1)%TAULEN] +
																FIELDhalf[(x*POSLEN+y)*TAULEN + (t-1+TAULEN)%TAULEN]) +

											  (1./(DX*DX))*(-2.*FIELDhalf[(x*POSLEN+y)*TAULEN + t] +                      //finite difference in x
											  					FIELDhalf[(x*POSLEN+(y+1)%POSLEN)*TAULEN + t] +
																FIELDhalf[(x*POSLEN+(y-1+POSLEN)%POSLEN)*TAULEN + t]) +

											  (1./(DX*DX))*(-2.*FIELDhalf[(x*POSLEN+y)*TAULEN + t] +                      //finite difference in y
											  					FIELDhalf[(((x-1+POSLEN)%POSLEN)*POSLEN+y)*TAULEN + t] +
																FIELDhalf[(((x+1)%POSLEN)*POSLEN+y)*TAULEN + t]);

					double HalfDrift = MASS*MASS*FIELDhalf[(x*POSLEN+y)*TAULEN + t] - PhiDoubPrimeHalf + INTERACTION*quadtrm*4*pow(FIELDhalf[(x*POSLEN+y)*TAULEN + t],3.);
					FIELDnext[(x*POSLEN+y)*TAULEN + t] = FIELDcurr[(x*POSLEN+y)*TAULEN + t] - dt5*0.5*(Drift[id]+HalfDrift) + sqrt(dt5)*NOISE[id];
				}

		std::swap(FIELDcurr,FIELDnext);

		//TEST SIN UPDATE
		for(int64_t x=0;x<POSLEN;x++)
  			for(int64_t y=0;y<POSLEN;y++)
				for(int64_t t=0;t<TAULEN;t++)
					FIELDcurr[(x*POSLEN+y)*TAULEN + t]=sin(2*PI*x/POSLEN)*sin(2*PI*y/POSLEN)*sin(2*PI*t/TAULEN);


		//WORKS <=> OMGLEN=TAULEN, MOMLEN=POSLEN
		for(int64_t i=0;i<OMGLEN*MOMLEN*MOMLEN;i++){
			FIELDOMGA[2*i]=FIELDcurr[i];
			FIELDOMGA[2*i+1]=0.;
		}

		fftw_execute(FTForInit);

		for(int64_t i=0;i<2*OMGLEN*MOMLEN*MOMLEN;i++)
			FIELDOMGA[i]/=sqrt(OMGLEN*MOMLEN*MOMLEN);

		//----------------------
		//BEGINNING MEASUREMENTS
		//----------------------
		if(t5>STARTMEAS){
			if(t5%MEASEACH==0){

				h++;
				//-------------------
				//EUCIDEAN CORRELATOR
				//-------------------
				for(int64_t i=0;i<TAULEN;i++)
					EuclCORRZero[i]=0.;

				for(int64_t x=0;x<POSLEN;x++)
					for(int64_t y=0;y<POSLEN;y++)
						for(int64_t t=0;t<TAULEN;t++){

							EuclCORR[(x*POSLEN+y)*TAULEN + t]=FIELDcurr[0]*FIELDcurr[(x*POSLEN+y)*TAULEN + t];
							EuclCORRZero[t]+=EuclCORR[(x*POSLEN+y)*TAULEN + t]/(POSLEN*POSLEN);
						}

				for(int64_t i=0;i<TAULEN;i++)
			  		EuclCORRZeroAvg[i]+=EuclCORRZero[i];

				//------------------
				//IM-FREQ CORRELATOR
				//------------------
				for(int64_t k=0;k<OMGLEN*BINLEN;k++){
	                ReCORRmom[k]=0.;
	                ImCORRmom[k]=0.;
	            }

				for(int64_t x=0,negx=0, PCNT=0; x<MOMLEN; x++, negx=MOMLEN*OMGLEN*(MOMLEN-x))
					for (int64_t y=0, negy=negx; y<MOMLEN; y++, negy=negx+OMGLEN*(MOMLEN-y))
						for (int64_t t=0, NCNT=negy; t<OMGLEN; t++, NCNT=negy+(OMGLEN-t),PCNT++){

							bindex = (int)floor(sqrt(momxsquared[x*MOMLEN+y])/DBIN);

							ReCORR[PCNT] = FIELDOMGA[2*PCNT]*FIELDOMGA[2*NCNT] - FIELDOMGA[2*PCNT+1]*FIELDOMGA[2*NCNT+1];
							ImCORR[PCNT] = FIELDOMGA[2*PCNT]*FIELDOMGA[2*NCNT+1] + FIELDOMGA[2*NCNT]*FIELDOMGA[2*PCNT+1];

							if(bindex<BINLEN){
								ReCORRmom[bindex*OMGLEN+t] += ReCORR[PCNT];
								ImCORRmom[bindex*OMGLEN+t] += ImCORR[PCNT];
							}
						}

				for(int64_t y=0;y<BINLEN;y++)
					for(int64_t x=0;x<OMGLEN;x++){
						if(mombin[y]!=0){
							ReCORRmom[y*OMGLEN+x]/=mombin[y];
							ImCORRmom[y*OMGLEN+x]/=mombin[y];
						}
						ReCORRmomAvg[y*OMGLEN+x]+=ReCORRmom[y*OMGLEN+x];
						ImCORRmomAvg[y*OMGLEN+x]+=ImCORRmom[y*OMGLEN+x];
					}

			}//almost end of measurements
		}//end of measurements

		//------------------------------------
		//OUTPUT FOR BAYESIAN ANALYSIS REPRISE
		//------------------------------------
		if( MEASBAYESCORREACH!=0 && (t5>STARTMEAS&&t5%MEASBAYESCORREACH==0) ){

			CORROUTI.write ((char*) ReCORR, sizeof(double)*OMGLEN);
			CORROUTE.write ((char*) EuclCORRZero, sizeof(double)*TAULEN);
			if(t5%(MEASBAYESCORREACH*10)==0){
				CORROUTI.flush();
				CORROUTE.flush();
			}
		}
	}//end of langevin

	CORROUTE.close();
	CORROUTI.close();

	//PRINT MEAN PROPAGATORS
    FILE* tizio;
    sprintf(fname,"./risultati/IMFREQmean");
    tizio = fopen(fname,"w");
    for(int64_t k=0;k<OMGLEN;k++)
        fprintf(tizio,"%ld %f %e %e\n",k,momtem[k],ReCORRmomAvg[k]/h,ImCORRmomAvg[k]/h);
    fclose(tizio);

	FILE* caio;
    sprintf(fname,"./risultati/EUCLmean");
    caio = fopen(fname,"w");
    for(int64_t t=0;t<TAULEN;t++)
    	fprintf(caio,"%d %e\n",t,EuclCORRZeroAvg[t]/h);
    fclose(caio);


}//end of main
