//--------------------------------//
//BAYESIAN SPECTRAL RECONSTRUCTION//
//--------------------------------//
//PRIOR MODEL m(I)=1

#include <fftw3.h>
#include <gmp.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <mpfr.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define PI (3.141592653589793)
// double Alpha = 0.5;

struct MinAux{

  double* AvgData;
  double* DiagCorr;
  double* Basis;
  double* Kernel;
  int Nomega;
  int Nmu;
  double Dmu;
  double Prior;
  double Alpha;
};

void PrintMatrix(double *A, int64_t Ncol){

    for(int i=0;i<5;i++){
        for(int j=0;j<5;j++)
            printf("%le ",A[i*Ncol+j]);
        printf("\n");
    }
}

void PrintGSLMatrix(gsl_matrix *A){

    for(int i=0;i<5;i++){
        for(int j=0;j<5;j++)
            printf("%lf ",gsl_matrix_get(A,i,j  ));
        printf("\n");
    }
}

//--------------------------------
//MAXIMIZING POSTERIOR PROBABILITY
//--------------------------------
double PostProbVals(const gsl_vector * x, void * params){			//f

    struct MinAux *Auxiliaries;
    Auxiliaries=(struct MinAux *)params;

    double *DiagCorr;
    double *Basis;
    double *AvgData;
    double *Kernel;
    int Nomega;
    int Nmu;
    double Dmu;
	double Prior;
	double Alpha;

    DiagCorr=(*Auxiliaries).DiagCorr;
    Basis=(*Auxiliaries).Basis;
    AvgData=(*Auxiliaries).AvgData;
    Kernel=(*Auxiliaries).Kernel;
    Nomega=(*Auxiliaries).Nomega;
    Nmu=(*Auxiliaries).Nmu;
    Dmu=(*Auxiliaries).Dmu;
	Prior=(*Auxiliaries).Prior;
	Alpha=(*Auxiliaries).Alpha;

	double *ModelData = calloc(Nomega,sizeof(double));

    double PostProb = 0;

    //CALCULATE MODEL DATA
    for(int64_t w=0;w<Nomega;w++)
        for(int64_t m=0;m<Nmu;m++)
            ModelData[w] += Dmu*Kernel[w*Nmu+m]*exp(gsl_vector_get(x,m))*Prior;

    //------------------------------------------------
    //FIRST PART OF POSTERIOR PROBABILITY - LIKELIHOOD
    //------------------------------------------------
    double Likelihood = 0;

    for(int64_t w=0;w<Nomega;w++)
        Likelihood += 0.5*pow(AvgData[w] - ModelData[w],2)/DiagCorr[w];

    //------------------------------------------------
    //SECOND PART OF POSTERIOR PROBABILITY - REGULATOR
    //------------------------------------------------
    double Regulator = 0;

    for(int64_t m=0;m<Nmu;m++)
        Regulator += Alpha*Dmu*(1 - Prior*exp(gsl_vector_get(x,m))/Prior + gsl_vector_get(x,m));

    PostProb = Likelihood - Regulator;

    free(ModelData);

	return PostProb;

}

void PostProbDeriv(const gsl_vector * x, void * params, gsl_vector * g){			//df

    struct MinAux *Auxiliaries;
    Auxiliaries=(struct MinAux *)params;

    double *DiagCorr;
    double *Basis;
    double *AvgData;
    double *Kernel;
    int Nomega;
    int Nmu;
    double Dmu;
	double Prior;
	double Alpha;

    DiagCorr=(*Auxiliaries).DiagCorr;
    Basis=(*Auxiliaries).Basis;
    AvgData=(*Auxiliaries).AvgData;
    Kernel=(*Auxiliaries).Kernel;
    Nomega=(*Auxiliaries).Nomega;
    Nmu=(*Auxiliaries).Nmu;
    Dmu=(*Auxiliaries).Dmu;
	Prior=(*Auxiliaries).Prior;
	Alpha=(*Auxiliaries).Alpha;

    double *ModelData = calloc(Nomega,sizeof(double));

    //CALCULATE MODEL DATA
    for(int64_t w=0;w<Nomega;w++)
        for(int64_t m=0;m<Nmu;m++)
            ModelData[w] += Dmu*Kernel[w*Nmu+m]*exp(gsl_vector_get(x,m))*Prior;

    for(int64_t k=0;k<Nmu;k++){

        double dPostProb = 0;

        //-------------------------------
        //CALCULATE LIKELIHOOD DERIVATIVE
        //-------------------------------
        double dLikelihood = 0;

        for(int64_t w=0;w<Nomega;w++)
            dLikelihood += -1*(AvgData[w] - ModelData[w])*(Dmu*Kernel[w*Nmu+k]*exp(gsl_vector_get(x,k))*Prior)/DiagCorr[w];

        //------------------------------
        //CALCULATE REGULATOR DERIVATIVE
        //------------------------------
        double dRegulator = 0;

//        for(int64_t m=0;m<Nmu;m++)
        dRegulator = Alpha*Dmu*(1 - Prior*exp(gsl_vector_get(x,k)));

        dPostProb = dLikelihood - dRegulator;
        gsl_vector_set(g,k,dPostProb);
    }

    free(ModelData);
}

void PostProbValsAndDeriv(const gsl_vector * x, void * params, double * f, gsl_vector * g){	//fdf

    *f = PostProbVals(x,params);
    PostProbDeriv(x,params,g);
}

//*******************//
//*******************//
//***MAIN*FUNCTION***//
//*******************//
//*******************//


int main(int argc, char *argv[]){

    //READING INITIAL CONDITIONS
    FILE* input, *output;
    //double k = atof(argv[1]);
    int64_t Nomega = 64;
	int64_t Nmu = 1000;
	int64_t Nmeas = 2000;
	double muMax = 2.2;
    double K = 1e-3;
    //char init[256];

    double Dmu = muMax/Nmu;
	double Prior = 1.;
	double Alpha = 0.;

	double* SpecFct = malloc(sizeof(double)*Nmu);
	double* Kernel = malloc(sizeof(double)*Nmu*Nomega);
    int* Tau = malloc(sizeof(int)*Nomega);

    double* AvgCorr = malloc(sizeof(double)*Nomega);
    double* RawCorr = malloc(sizeof(double)*Nomega*Nmeas);
    double* ModelData = calloc(Nomega,sizeof(double));
    double* momtem = malloc(sizeof(double)*Nomega);
    double* Cov_Matrix = calloc(Nomega*Nomega,sizeof(double));
    double* Cov_Basis = malloc(sizeof(double)*Nomega*Nomega);
    double* Cov_Eigenval = malloc(sizeof(double)*Nomega);

    char* fname = malloc(sizeof(char)*256);
    char dirname[256];

    sprintf(dirname,"risultati-imfreq-free-64");

    //-----------------------------
    //SETUP RANDOM NUMBER GENERATOR
    //-----------------------------
    const gsl_rng_type * A;
    A = gsl_rng_default;
    gsl_rng_env_setup();
    gsl_rng *rnd = gsl_rng_alloc(A);

    //------------
    //READING DATA
    //------------

// 	sprintf(fname,"/remote/pi326b/rizzardi/symm-imfreq-free/symm-imfreq-freeIMFREQmean");
// 	input = fopen(fname,"r");
// 	for(int i=0;i<Nomega;i++){
// 		int a;
// 		double b;
//         fscanf(input,"%d %lf %lf %lf\n",&a,&momtem[i],&AvgCorr[i],&b);
//     }
// 	fclose(input);


    sprintf(fname,"./momtem");
    input = fopen(fname,"r");
    for(int i=0;i<Nomega;i++){
        fscanf(input, "%lf\n",&momtem[i]);
    }
    fclose(input);

    sprintf(fname,"./aaa.txt");
    input = fopen(fname,"r");
    for(int i=0;i<Nomega;i++){
        fscanf(input, "%lf\n",&AvgCorr[i]);
    }
    fclose(input);


//     for(int64_t t=0;t<Nmeas;t++){
//         int a;
//         double b,c;
// 		printf("Reading file %d/%d\r",t,Nmeas);
//         sprintf(fname,"/remote/pi326b/rizzardi/symm-imfreq-free/symm-imfreq-freeIMFREQmean.%d",t+39900);
//         input = fopen(fname,"r");
// 		for(int64_t w=0;w<Nomega;w++)
// 			fscanf(input,"%d %lf %lf %lf\n",&a,&b,&RawCorr[w*Nmeas+t],&c);
//         fclose(input);
//     }




    //ADDING GAUSSIAN NOISE

    for(int64_t j=0;j<Nomega;j++){
        for(int64_t k=0;k<Nmeas;k++){
            double sigma = AvgCorr[j]*K*sqrt(Nmeas);
            RawCorr[j*Nmeas + k] = AvgCorr[j] + gsl_ran_gaussian(rnd, sigma);///(2*PI*sigma*sigma);
        }
    }

    //--------------------------------------------------------------------------------
    //AUTOCORRELATION COMPUTATION - ALREADY DONE IN PYTHON - just reading coefficients
    //--------------------------------s------------------------------------------------
    printf("Computing autocorrelation... ");

    fftw_complex *obs = fftw_malloc( sizeof(fftw_complex)*Nmeas);

    fftw_plan FFT, IFFT;
    FFT= fftw_plan_dft_1d(Nmeas, obs,obs,FFTW_FORWARD,FFTW_ESTIMATE);
    IFFT= fftw_plan_dft_1d(Nmeas, obs,obs,FFTW_BACKWARD,FFTW_ESTIMATE);

    for(int64_t w=0;w<Nomega;w++){
        for(int64_t k=0;k<Nmeas;k++){

            obs[k][0]=RawCorr[w*Nmeas+k] - AvgCorr[w];
            obs[k][1]=0.;
        }

        // CARRY OUT THE COMPUTATION OF THE AUTOCORRELATION VIA FOURIER TRANSFORM

        fftw_execute(FFT);

        for(int k=0;k<Nmeas;k++) {
            double abssqr= obs[k][0]*obs[k][0] + obs[k][1]*obs[k][1];
            obs[k][0]=abssqr/Nmeas;
            obs[k][1]=0.;
        }

        fftw_execute(IFFT);

        for(int64_t t=1;t<Nmeas;t++){
            double autoc = obs[t][0]/obs[0][0];
            if(autoc < 0.08){
                Tau[w] = t;
                break;
            }
        }
    }
    printf("OK\n");

    //----------------------------------------------------------------
    //COMPUTING COVARIANCE MATRIX - same notation as Jarrel-Gubernatis
    //----------------------------------------------------------------
    printf("Computing covariance matrix...");
    for(int64_t i=0;i<Nomega;i++)
        for(int64_t k=0;k<Nomega;k++)
            for(int64_t j=0;j<Nmeas;j++)
                Cov_Matrix[i*Nomega+k] += Tau[i]*Tau[k]*(AvgCorr[i]-RawCorr[i*Nmeas+j])*(AvgCorr[k]-RawCorr[k*Nmeas+j])/(Nmeas*(Nmeas-1));
    printf(" OK\n");

	PrintMatrix(Cov_Matrix,Nomega);
    
	//
    //---------------------------------
    //DIAGONALIZATION COVARIANCE MATRIX
    //---------------------------------
    printf("Diagonalization covariance matrix...");

    //need to use gsl_matrix and gsl_vector. simple casting does not work
    gsl_eigen_symmv_workspace* Eigensystem = gsl_eigen_symmv_alloc (Nomega);

    gsl_matrix* CovMatrix = gsl_matrix_calloc(Nomega,Nomega);
    gsl_matrix* CovBasis = gsl_matrix_calloc(Nomega,Nomega);
    gsl_vector* CovEigenval = gsl_vector_alloc(Nomega);
    gsl_vector* Avg_Corr = gsl_vector_alloc(Nomega);
    gsl_vector* Mult = gsl_vector_alloc(Nomega);

    for(int64_t i=0;i<Nomega;i++)
        for(int64_t j=0;j<Nomega;j++)
            gsl_matrix_set(CovMatrix,i,j,Cov_Matrix[i*Nomega+j]);

    for(int64_t i=0;i<Nomega;i++)
        gsl_vector_set(Avg_Corr,i,AvgCorr[i]);

    gsl_eigen_symmv(CovMatrix, CovEigenval, CovBasis, Eigensystem);

    for(int64_t i=0;i<Nomega;i++)
        Cov_Eigenval[i] = gsl_vector_get(CovEigenval,i);

    gsl_blas_dgemv(CblasNoTrans, 1, CovBasis, Avg_Corr, 0, Mult);

    for(int64_t i=0;i<Nomega;i++)
        AvgCorr[i] = gsl_vector_get(Mult,i);

    for(int64_t i=0;i<Nomega;i++)
        for(int64_t k=0;k<Nomega;k++)
            Cov_Basis[i*Nomega+k] = gsl_matrix_get(CovBasis,i,k);

    gsl_eigen_symmv_free(Eigensystem);
    gsl_vector_free(CovEigenval);
    gsl_vector_free(Avg_Corr);
    gsl_vector_free(Mult);
    gsl_matrix_free(CovMatrix);
    free(Cov_Matrix);
    printf(" OK\n");

// 	 for(int64_t i=0;i<Nomega;i++)
//         printf("%lf\n",Cov_Eigenval[i]);

// 	 return 0;

    //------------------
    //KERNEL COMPUTATION
    //------------------
    printf("Computing kernel...\n");
    for(int64_t k=0;k<Nomega;k++) momtem[k]=sqrt(2.-2.*cos(k*2*PI/Nomega));
    for(int64_t w=0;w<Nomega;w++)
        for(int64_t m=0;m<Nmu;m++)
            Kernel[w*Nmu+m] = 2*m*Dmu*m*Dmu/(m*Dmu*m*Dmu+momtem[w]*momtem[w]);

    Kernel[0] = 0.;

    gsl_matrix* UncKernel = gsl_matrix_alloc(Nomega,Nmu);
    gsl_matrix* Result = gsl_matrix_calloc(Nomega,Nmu);

    for(int w=0;w<Nomega;w++)
        for(int m=0;m<Nmu;m++)
            gsl_matrix_set(UncKernel,w,m,Kernel[w*Nmu+m]);

    //COMPUTING "UNCORRELATED" KERNEL
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, CovBasis, UncKernel, 0, Result);

    for(int w=0;w<Nomega;w++)
        for(int m=0;m<Nmu;m++)
            Kernel[w*Nmu+m] = gsl_matrix_get(Result,w,m);

    gsl_matrix_free(CovBasis);
    printf("OK\n");

    //printf("%lf\n",Kernel[Nmu+1]);

    // output = fopen("kernel","w");
    // for(int w=0;w<Nomega;w++){
    //     for(int m=0;m<Nmu;m++)
    //         fprintf(output,"%lf ",Kernel[w*Nmu+m]);
    //     fprintf(output,"\n");
    // }
    // fclose(output);



    //-------------------
    //SETUP THE MINIMIZER
    //-------------------
    printf("Setting up the minimizer...\t");
    struct MinAux MinimizerAuxiliaries;

    MinimizerAuxiliaries.AvgData=AvgCorr;
    MinimizerAuxiliaries.DiagCorr=Cov_Eigenval;
    MinimizerAuxiliaries.Basis=Cov_Basis;
    MinimizerAuxiliaries.Kernel=Kernel;
    MinimizerAuxiliaries.Nomega=Nomega;
    MinimizerAuxiliaries.Nmu=Nmu;
    MinimizerAuxiliaries.Dmu=Dmu;
	MinimizerAuxiliaries.Prior=Prior;
	MinimizerAuxiliaries.Alpha=Alpha;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    T = gsl_multimin_fdfminimizer_vector_bfgs2;           //type of the minimizer
    s = gsl_multimin_fdfminimizer_alloc (T, Nmu);         //instance of the minimizer

    gsl_multimin_function_fdf Posterior_Prob;             //function to minimize

    gsl_vector *ExpFitParams = gsl_vector_calloc(Nmu);    //fit parameters

    Posterior_Prob.n = Nmu;
    Posterior_Prob.f = PostProbVals;
    Posterior_Prob.df = PostProbDeriv;
    Posterior_Prob.fdf = PostProbValsAndDeriv;
    Posterior_Prob.params = (void*) &MinimizerAuxiliaries;
    printf("OK\n");

    //--------------------
    //MINIMIZATION PROCESS
    //--------------------

	//find alpha with bisection
	int check = 0;
    int status = 0;
    double a=atof(argv[1]), b=atof(argv[2]),c;
    double f_a, f_b;
    double Qold, Qnew, RelChange;
	int64_t loop = 0;

    printf("Proceeding with minimization...\n");

    //----------------------------------------------------------//
    //BISECTION METHOD SOLVING Q + alpha*S - Nomega = 0 IN alpha//
    //EVERY FUNCTION EVALUATION REQUIRES MINIMIZATION-----------//
    //----------------------------------------------------------//

    //LOOP TO CHECK SOLVING CONDITION
    do
    {
        c = (a+b)/2;
        MinimizerAuxiliaries.Alpha = c;

        int step = 0;
		//INITIAL CONDITIONS IN PARAMETER SPACE
		for(int k=0; k<Nmu;k++)
			gsl_vector_set(ExpFitParams,k,1.);
		//SETUP OF THE MINIMIZER
		gsl_multimin_fdfminimizer_set (s, &Posterior_Prob, ExpFitParams, 0.01, 0.01);

        printf("Interval: [%f,%f] ",a,b);

		do
		{
            Qold = s->f;
            //status = GSL_CONTINUE;
            status = gsl_multimin_fdfminimizer_iterate (s);
            Qnew = s->f;
            RelChange = (Qold-Qnew)/Qold;
            if(step%50==0){
// 				double S=0.;
// 				for(int64_t m=0;m<Nmu;m++)
// 					S += Dmu*(1 - exp(gsl_vector_get(s->x,m)) + gsl_vector_get(s->x,m));
//                 double L = Qnew + (MinimizerAuxiliaries.Alpha)*S;
                printf("%le -- %le\n",gsl_blas_dnrm2(s->gradient),RelChange);
			}
            status = GSL_CONTINUE;
            if(RelChange<5e-9 && RelChange>0)
                status = GSL_SUCCESS;
            //status = GSL_CONTINUE;
            step++;

		}
		while (status == GSL_CONTINUE && step < 100000);

		//FINDING alpha
		double Q = s->f;
		double S=0;

		for(int64_t m=0;m<Nmu;m++)
			S += Dmu*(1 - exp(gsl_vector_get(s->x,m)) + gsl_vector_get(s->x,m));

		double L = Q + (MinimizerAuxiliaries.Alpha)*S;

        printf("- c = %le - |f(c)|: %e - L:  %f\n",c,fabs(L-Nomega),L);

        if(L-Nomega<0)
            a=c;
        else
            b=c;


        if(fabs(L-Nomega) < 5e-1)
            check = 1;

    }//end condition loop
    while (check == 0);

    printf(" OK\n");
    printf("\nAlpha: %f\n",c);

    printf("Computing spectral function...");
    for(int64_t i=0;i<Nmu;i++)
        SpecFct[i] = i*Dmu*Prior*exp(gsl_vector_get(s->x,i));

    sprintf(fname,"SpecFct");
	output = fopen(fname,"w");
	for(int64_t i=0;i<Nmu;i++)
		fprintf(output,"%f %e\n",i*Dmu,SpecFct[i]);
	fclose(output);
    printf(" Success!\n");


    fftw_free(FFT);
    fftw_free(IFFT);

    free(ModelData);
    free(RawCorr);
    free(AvgCorr);

    gsl_multimin_fdfminimizer_free(s);

	return 0;
}//end main
