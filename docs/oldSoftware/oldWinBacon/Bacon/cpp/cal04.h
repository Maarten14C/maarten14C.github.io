


#ifndef CAL_H
#define CAL_H 

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include "ranfun.h"
#include "Matrix.h"

#define IntCal04ROWS 3309
#define IntCal04COLS 3


class Cal {

protected:

	int k;
	double mu, sig;

public: 
	
	Cal(int kk) { k = kk; }
	
	double GetSig() { return sig; }
	double GetMu() { return mu; }
	
	virtual char *Name() = 0;
	
	virtual double cal(double theta) = 0;               
	virtual double U( double y, double vr, double theta) = 0;    
	virtual double Ut( double y, double vr, double theta, double a, double b) = 0;    
	
	virtual double MinCal() = 0;
	virtual double MaxCal() = 0;

};


//Constant "calibration curve", no curve, basically
class ConstCal : public Cal {

public:

	ConstCal() : Cal(0) {
		printf("Constant calibration curve.\n");
	}
	
	double cal(double theta) {
		mu = theta;
		sig = 0.0;
		
		return theta;
	}

	char *Name() { return "Constant c. curve"; }

	double U( double y, double vr, double theta) { return 0.5*sqr(y - theta)/vr; }    
	double Ut( double y, double vr, double theta, double a, double b) { return (a + 0.5)*log( b + 0.5*sqr(y - theta)/vr); }

	double MinCal() { return -1.0e300; }
	double MaxCal() { return 1.0e300; }
};



/**** IntCal04	WITH THE BOMB PEAK ADDED, INTERPOLATED EVERY 5 YEARS FROM -45 TO -10 ****/

/*
 IntCal09 Northern Hemisphere atmospheric radiocarbon calibration curve
 Format: cal BP, 14C BP, ± Error is age-corrected as per Stuiver and Polach (1977).
 The model used for curve construction is presented in Heaton et al. (2009).
 Note that the data spacing changes from 5 years for the range from 0 to 11.2 to cal kBP, 10 yrs
 for 11.2–15 cal kBP, 20 yrs for 15–25 cal kBP, 50 yrs for 25–40 cal kBP, and 100 yrs for 40–50 cal kBP.  
*/
class IntCal04 : public Cal {

protected:

	Matrix *CCB;
	SubMatrix CC;
	
	SubMatrix A;
	
	double const2;
	
public:

	IntCal04(char *fnam) : Cal(IntCal04ROWS) {
	
		CCB = new Matrix( IntCal04ROWS, IntCal04COLS);
		
		CC.Set( CCB, CCB->nRow(), CCB->nCol());
		
		printf("IntCal04, with bomb peak added: Reading from file: %s\n", fnam);
		
		if (CC.filescan(fnam) == 0) {
			printf("Cal: ERROR: Could not find IntCal04 cal. curve, file not found: %s\n", fnam);
			
			exit(0);
		}
		
		const2 = 0.5*log(2.0 * M_PI); //M_PI defined in gsl library
	
		/*A.Set( CCB, 20, CCB->nCol());
				
		A.print("IntCal04:\n");
		printf("...\n");*/
		
	}
	
	~IntCal04() {
		A.~SubMatrix();
		CC.~SubMatrix();
		delete CCB;
	}
	

	double MinCal() { return -45.0; }
	double MaxCal() { return 25980.0; }

	virtual char *Name() { return "IntCal04"; }
	
	
//prototype for gsl_compare
// x==y, 0, x<y, -1, x>y, 1.
//int fcmp (double x, double y, double epsilon = 0.00000000001); 
	double cal(double theta)
	{
        if (fcmp(theta, -45.0) == -1) //******NB: -45, starting point of current curve      
        {
                //fprintf( stderr, "WARNING: Calibration attempted beyond IntCal04 cal. curve limits, theta= %f\n",theta);
                k = 0;
                mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5;
                sig <- CC(k,2);
        }
        else
                if (fcmp(theta, 12500.0) != 1)
                {				//************** NB: 9 is the number of rows before cal year 0 (in offical IntaCal04 this is 1, and so forth)
                        k = 1 + 8 + (int) floor(theta/5.0);
                        mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5.0;
                        sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/5.0;
                }
                else
                        if (fcmp(theta, 15000.0) != 1)
                        {
                                k = 2501 + 8 + (int) floor((theta-12500.0)/10.0);
                                mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/10.0;
                                sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/10.0;
                        }
                        else
                                if (fcmp(theta, 25980.0) != 1)
                                {
                                        k = 2751 + 8 + (int) floor((theta-15000.0)/20.0);
                                        mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/20.0;
                                        sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/20.0;
                                }
                                else
                                {
                        //fprintf( stderr, "WARNING: Calibration attempted beyond IntCal04 cal. curve limits, theta= %f\n",theta);
                                        k = 3299 + 8;
                                        mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/20.0;
                                        sig =CC(k,2);
                                }
        
		return mu;
	}
	
	
	               
	double U( double y, double vr, double theta)                
	{
        cal(theta);
        
        double tau = 1.0/(vr + sqr(sig));
        
        return const2 - 0.5*log(tau)  + tau * 0.5 * sqr(y-mu);
	}

	
	//The new energy using the t distribution
	double Ut( double y, double vr, double theta, double a, double b)                
	{
        cal(theta);
        
        double tau = 1.0/(vr + sqr(sig));
        
        return (a + 0.5)*log( b + tau * 0.5 * sqr(y-mu));
	}


};



//The Marine curve of IntCal04 has the same format as the terrestrial one.	

/*
 Note that the data spacing changes from 5 years for the range from 0 to 12.5 to cal kBP,
 10 yrs for 11.2–15 cal kBP, 20 yrs for 15–25 cal kBP, 50 yrs for 25–40 cal kBP, and 100 yrs for 40–50 cal kBP
*/
class IntCal04Marine : public Cal {

protected:

	Matrix *CCB;
	SubMatrix CC;
	
	SubMatrix A;
	
	double const2;
	
public:

	IntCal04Marine(char *fnam) : Cal(3301) {
	
		CCB = new Matrix( 3301, 3);
		
		CC.Set( CCB, CCB->nRow(), CCB->nCol());
		
		printf("IntCal04Marine: Reading from file: %s\n", fnam);
		
		if (CC.filescan(fnam) == 0) {
			printf("Cal: ERROR: Could not find IntCal04Marine cal. curve, file not found: %s\n", fnam);
			
			exit(0);
		}
		
		const2 = 0.5*log(2.0 * M_PI); //M_PI defined in gsl library
	
		/*A.Set( CCB, 20, CCB->nCol());
				
		A.print("IntCal04Marine:\n");
		printf("...\n");*/
		
	}
	
	~IntCal04Marine() {
		A.~SubMatrix();
		CC.~SubMatrix();
		delete CCB;
	}
	

	double MinCal() { return 0.0; }
	double MaxCal() { return 25980.0; }


		char *Name() { return "IntCal04Marine"; }
	
	
//prototype for gsl_compare
// x==y, 0, x<y, -1, x>y, 1.
//int fcmp (double x, double y, double epsilon = 0.00000000001); 
	double cal(double theta)
	{
        if (fcmp(theta, 0.0) == -1) //******NB: -45, starting point of current curve      
        {
                //fprintf( stderr, "WARNING: Calibration attempted beyond IntCal04 cal. curve limits, theta= %f\n",theta);
                k = 0;
                mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5;
                sig <- CC(k,2);
        }
        else
                if (fcmp(theta, 12500.0) != 1)
                {				//************** NB: 0 is the number of rows before cal year 0 
                        k = 0 + (int) floor(theta/5.0);
                        mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5.0;
                        sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/5.0;
                }
                else
                        if (fcmp(theta, 15000.0) != 1)
                        {
                                k = 2500 + (int) floor((theta-12500.0)/10.0);
                                mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/10.0;
                                sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/10.0;
                        }
                        else
                                if (fcmp(theta, 26000.0) != 1)
                                {
                                        k = 2750 + (int) floor((theta-15000.0)/20.0);
                                        mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/20.0;
                                        sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/20.0;
                                }
                                else
                                {
                        //fprintf( stderr, "WARNING: Calibration attempted beyond IntCal04 cal. curve limits, theta= %f\n",theta);
                                        k = 3299;
                                        mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/20.0;
                                        sig =CC(k,2);
                                }
        
		return mu;
	}
	
	
	               
	double U( double y, double vr, double theta)                
	{
        cal(theta);
        
        double tau = 1.0/(vr + sqr(sig));
        
        return const2 - 0.5*log(tau)  + tau * 0.5 * sqr(y-mu);
	}

	
	//The new energy using the t distribution
	double Ut( double y, double vr, double theta, double a, double b)                
	{
        cal(theta);
        
        double tau = 1.0/(vr + sqr(sig));
        
        return (a + 0.5)*log( b + tau * 0.5 * sqr(y-mu));
	}

};




/**** Southern Hemisphere c curve ****/
class SHCal04 : public Cal {

protected:

	Matrix *CCB;
	SubMatrix CC;
	
	SubMatrix A;
	
	double const2;
	
public:

	SHCal04(char *fnam) : Cal(2482) {
	
		CCB = new Matrix( 2482, 3);
		
		CC.Set( CCB, CCB->nRow(), CCB->nCol());
		
		printf("SHCal04: Reading from file: %s\n", fnam);
		
		if (CC.filescan(fnam) == 0) {
			printf("Cal: ERROR: Could not find SHCal04 cal. curve, file not found: %s\n", fnam);
			
			exit(0);
		}
		
		const2 = 0.5*log(2.0 * M_PI); //M_PI defined in gsl library
	
		/*A.Set( CCB, 20, CCB->nCol());
				
		A.print("IntCal04:\n");
		printf("...\n");*/
		
	}
	
	~SHCal04() {
		A.~SubMatrix();
		CC.~SubMatrix();
		delete CCB;
	}
	

	double MinCal() { return -5.0; }
	double MaxCal() { return 12400.0; }

	virtual char *Name() { return "SHCal04"; }
	
	
//prototype for gsl_compare
// x==y, 0, x<y, -1, x>y, 1.
//int fcmp (double x, double y, double epsilon = 0.00000000001); 
	double cal(double theta)
	{
        if (fcmp(theta, -5) == -1) //******NB: -45, starting point of current curve      
        {
                //fprintf( stderr, "WARNING: Calibration attempted beyond IntCal04 cal. curve limits, theta= %f\n",theta);
                k = 0;
                mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5;
                sig <- CC(k,2);
        }
        else
                if (fcmp(theta, 12400.0) != 1)
                {				//************** NB: 1 is the number of rows before cal year 0 
                        k = 1 + (int) floor(theta/5.0);
                        mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5.0;
                        sig =CC(k,2) + (theta-CC(k,0))*(CC(k+1,2)-CC(k,2))/5.0;
                }
                else
				{
                        //fprintf( stderr, "WARNING: Calibration attempted beyond IntCal04 cal. curve limits, theta= %f\n",theta);
					k = 2480;
					mu = CC(k,1) + (theta-CC(k,0))*(CC(k+1,1)-CC(k,1))/5.0;
					sig =CC(k,2);
				}
        
		return mu;
	}
	
	
	               
	double U( double y, double vr, double theta)                
	{
        cal(theta);
        
        double tau = 1.0/(vr + sqr(sig));
        
        return const2 - 0.5*log(tau)  + tau * 0.5 * sqr(y-mu);
	}

	
	//The new energy using the t distribution
	double Ut( double y, double vr, double theta, double a, double b)                
	{
        cal(theta);
        
        double tau = 1.0/(vr + sqr(sig));
        
        return (a + 0.5)*log( b + tau * 0.5 * sqr(y-mu));
	}


};




/******* Class Det *******/


//Definition for class Det, to hold a generic "determination", could be non-radiocarbon
class Det {
protected:
	char *nm;	// lab number
	double y;		// mean
	double std;	// std dev

	double x;  //depth or anything else
		
	double deltaR;	// delta-R (reservoir correction)
	double deltaSTD;	// delta-R std dev. (reservoir correction)
	
	double a, b; //prior parameters for the t distribution

	Cal *cc;	// calibration curve to use
		
	double med;	// mean and variance and std dev. after reservoir correction
	double vr;
	double corrstd;

public:
	Det(char *enm, double ey, double estd, double xx, double edeltaR, double edeltaSTD, double ea, double eb, Cal *ecc) {

		//Read members
		nm = strdup(enm);
		y = ey;
		std = estd;
		x = xx;
		deltaR = edeltaR;
		deltaSTD = edeltaSTD;
		
		a = ea;
		b = eb;

		cc = ecc;

		med = y - deltaR;
		vr = sqr(std) + sqr(deltaSTD);
		corrstd = sqrt(vr);

	}
		
	void ShortOut() {
		printf("%s: %6.0f+-%-6.0f  d=%-g  ResCorr= %6.1f+-%-6.1f  a=%-g b=%-g   cc=%s\n",
				nm, y, std, x, deltaR, deltaSTD, a, b, cc->Name());
	}

	const char *labnm() { return nm; }
	double mean() { return y; }
	double sd() { return std; }
	double corr_mean() { return med; }
	double corr_vr() { return vr; }
	double res_mean() { return deltaR; }
	double res_std() { return deltaSTD; }
	double d() { return x; }

	//exp(-U) will be the likelihood for this determination.
	double U(double theta) { return cc->U( med, vr, theta); }
	double Ut(double theta) { return cc->Ut( med, vr, theta, a, b); }
};


//To hold a series of determinations
class Dets {

protected:

	Det **det;  //array with the determinations
	
	int m; //current number of determinations
	int max_m; //Maximum number of determinations

public:
	//The constructor only opens an array of pointers to Det structures
	Dets(int emax_m) {
	
		max_m = emax_m;
		
		m = 0;
		
		det = new Det * [max_m];
	}
	
	void AddDet(Det *de) {
	
		if (m == max_m) {
		
			printf("ERROR: Maximum number of determinations exceeded\n\n");
			
			exit(0);
		}

		m++;
		det[m-1] = de;

		printf("Added det: ");
		ShortOut(m-1);
	}
	
	int Size() { return m; }
	
	void ShortOut(int j) { det[j]->ShortOut(); } 
	
	const char *labnm(int j) { return det[j]->labnm(); }
	double y(int j)  {  return det[j]->corr_mean(); }
	double sd(int j) {  return det[j]->sd(); }
	double vr(int j) {  return det[j]->corr_vr(); }
	double d(int j)  {  return det[j]->d(); }
	

	double U(int j, double theta)  { return det[j]->U(theta); }
	double Ut(int j, double theta) { return det[j]->Ut(theta); }
};




#endif

