// This is the main DLL file.

#include "stdafx.h"
#include "stdio.h"
#include "math.h"

#include "IsotopomerDynamics.h"

using namespace System;
using namespace IsotopomerDynamics;

float pH = 0.00015574;  // pH natural abundance of 2H in nature;

int MassIsotopomers::CalculateMIDynamics(array <float, 1> ^ NaturalIsotopes, 
	array <float, 1> ^ fLabeledIsotopes, float fBWE1, int NEH, int Nall_Hydrogens)
{

	float a1, a2, a3;
	
	float pX_t, pX_H_ratio, fSpectShift;

	pX_t = fBWE1;

	pX_H_ratio = 1. - pX_t / (1. - pH);

	a1 = (float)NEH * (pX_t + pH) / (1. - pH - pX_t);

	a2 = (float)NEH * ((float)NEH - 1.) / 2.;

	a2 = a2 * Math::Pow((pX_t + pH) / (1. - pH - pX_t), 2.0);

	a3 = (float)NEH * ((float)NEH - 1.) *  ((float)NEH - 2.) / 6.;

	a3 = a3 * Math::Pow((pX_t + pH) / (1. - pH - pX_t), 3.);

	

	// the monoisotope
	fLabeledIsotopes[0] = Math::Pow(pX_H_ratio, NEH) * NaturalIsotopes[0];
	
	fLabeledIsotopes[1] = Math::Pow(pX_H_ratio, NEH - 1) * pX_t * NEH * NaturalIsotopes[0] / Math::Pow(1. - pH, 2);

	fLabeledIsotopes[1] = fLabeledIsotopes[1] + Math::Pow(pX_H_ratio, NEH) * NaturalIsotopes[1];

	fLabeledIsotopes[2] = NaturalIsotopes[2] / NaturalIsotopes[0] - (NaturalIsotopes[1] / NaturalIsotopes[0]) * pH * (float)NEH / (1. - pH);

	fLabeledIsotopes[2] = fLabeledIsotopes[2] + Math::Pow(pH / (1. - pH), 2.) * (float)NEH * ((float)NEH + 1.) / 2.0;

	fLabeledIsotopes[2] = fLabeledIsotopes[2] - Math::Pow((pX_t + pH) / (1 - pH - pX_t), 2) * (float)NEH * ((float)NEH + 1.) / 2.0;

	fLabeledIsotopes[2] = fLabeledIsotopes[2] * fLabeledIsotopes[0];

	fLabeledIsotopes[2] = fLabeledIsotopes[2] + (float)NEH * (pX_t + pH) * fLabeledIsotopes[1] / (1. - pH - pX_t);


	fLabeledIsotopes[3] = a3 * fLabeledIsotopes[0] + a2 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) +

				 a1 * (fLabeledIsotopes[2] - a1 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) - a2 * fLabeledIsotopes[0]);


	fLabeledIsotopes[3] = fLabeledIsotopes[3] + fLabeledIsotopes[0] * (NaturalIsotopes[3] / NaturalIsotopes[0] - (float)NEH * pH / (1. - pH) *
		(NaturalIsotopes[2] / NaturalIsotopes[0]  - ((float)NEH + 1.) * pH / 2. / (1. - pH) * 
		(NaturalIsotopes[1] / NaturalIsotopes[0] - (float)Nall_Hydrogens * pH / (1. - pH) ) - 
			(float) (3 * Nall_Hydrogens * NEH + 3 * Nall_Hydrogens - NEH * NEH - 3 * NEH - 2) *
			Math::Pow(pH / (1. - pH), 2)/6. ) );


	fLabeledIsotopes[4] = 1. - (fLabeledIsotopes[0] + fLabeledIsotopes[1] + fLabeledIsotopes[2] +
		                fLabeledIsotopes[3]);

	// computations for I4
	// I4 is a sum of five terms: p(Rest = 4) * p (X_h = 0) + p(Rest = 3) * p (X_H = 1) + ... + p(Rest = 0) * p(X_H = 4)

	//computing p(Rest = 0) * p(X_H = 4)
	float f_rest0_XH_4 = 0.;

	f_rest0_XH_4 = (float)(NEH * (NEH - 1)* (NEH - 2) * (NEH - 3) / 4 / 3 / 2) * powf((pH + pX_t) / (1. - pH - pX_t), 4.);

	f_rest0_XH_4 = f_rest0_XH_4 * fLabeledIsotopes[0];

	//end of computing p(Rest = 0) * p(X_H = 4)

	//computing p(Rest = 1) * p(X_H = 3). It is computed using:
	// p(Rest = 1) * p(X_H = 3) = p(Rest = 1) * p(X_H = 0) * p(X_H = 3)/ p(X_H = 0)
	// this will use a term from I1(t);

	float f_rest1_XH_3 = 0.;
	//p(Rest = 1) * p(X_H = 0) is:

	f_rest1_XH_3 = powf(pX_H_ratio, NEH) * NaturalIsotopes[1] -

		powf(pX_H_ratio, NEH) * NEH * pH / (1. - pH) * NaturalIsotopes[0];

	// multiply p(Rest = 1) * p(X_H = 0) to p(X_H = 3)/ p(X_H = 0)

	f_rest1_XH_3 = f_rest1_XH_3 * powf((pX_t + pH) / (1. - pX_t - pH), 3);

	f_rest1_XH_3 = f_rest1_XH_3 * (float)(NEH * (NEH - 1)*(NEH - 2)) / 3. / 2.;


	// end of computing p(Rest = 1) * p(X_H = 3)


	// computing p(Rest = 2) * p(X_H = 2). 
	// It is computed using:
	// p(Rest = 2) * p(X_H = 2) = p(Rest = 2) * p(X_H = 0) * p(X_H = 2)/ p(X_H = 0)
	// this will use a term from I2(t) for  p(Rest = 2) * p(X_H = 0);

	float f_rest2_XH_2 = 0.;

	//p(Rest = 2) * p(X_H = 0) is:

	f_rest2_XH_2 = NaturalIsotopes[2] - pH * NEH * NaturalIsotopes[1] / (1. - pH) +

		pH *pH * NEH * Nall_Hydrogens * NaturalIsotopes[0] / powf(1. - pH, 2.);

	f_rest2_XH_2 = f_rest2_XH_2 - (float) (NEH * (NEH + Nall_Hydrogens - 1)) * pH * pH *

		NaturalIsotopes[0] / powf(1. - pH, 2.) / 2.;

	f_rest2_XH_2 = f_rest2_XH_2 * powf(pX_H_ratio, NEH);

	// now multiply to p(X_H = 2)/ p(X_H = 0)

	f_rest2_XH_2 = f_rest2_XH_2 * powf((pX_t + pH) / (1. - pX_t - pH), 2);

	f_rest2_XH_2 = f_rest2_XH_2 *  (float)(NEH * (NEH - 1)) / 2.;

	// end of computing p(Rest = 2) * p(X_H = 2)


	// computing p(Rest = 3) * p(X_H = 1). 
	// It is computed using:
	// p(Rest = 3) * p(X_H = 1) = p(Rest = 3) * p(X_H = 0) * p(X_H = 1)/ p(X_H = 0)
	// this will use a term from I3(t) for  p(Rest = 3) * p(X_H = 1);
	float f_rest3_XH_1 = 0., f_rest3_XH_0 = .0;

	f_rest3_XH_0 = NaturalIsotopes[3] - (float) NEH * pH / (1. - pH) * (NaturalIsotopes[2] - 
		
		0.5 * (float) (NEH + 1) * pH / (1. - pH) * (NaturalIsotopes[1] - Nall_Hydrogens * pH *  NaturalIsotopes[0] / (1. - pH) ) -
		
		powf(pH / ( 1. - pH), 2.0) / 6. * (float) (Nall_Hydrogens * (NEH + 1) + (Nall_Hydrogens - 1) * (NEH + 2) +
		
			NEH * (Nall_Hydrogens - NEH - 2)) * NaturalIsotopes[0] );

	f_rest3_XH_0 = f_rest3_XH_0 * powf(pX_H_ratio, NEH);

	// now multiply to p(X_H = 1)/ p(X_H = 0) to obtain f_rest3_XH_1:

	f_rest3_XH_1 = f_rest3_XH_0 * (pX_t + pH) / (1. - pX_t - pH) * (float)NEH;

	// end of computing p(Rest = 3) * p(X_H = 1)

	
	// computing p(Rest = 4) * p(X_H = 0). 
	// It is computed using:
	// p(Rest = 4) = p(Rest-H = 4) * p(H=0) + p(Rest-H = 3) * p(H=1) + p(Rest-H=2)*p(H=2) +
	//               p(Rest-H = 1) * p(H = 3) + p(Rest-H = 0) * p(H = 4);
	//
	float f_rest4_XH_0 = 0;

	f_rest4_XH_0 = NaturalIsotopes[4] - (float) (NEH) * pH / (1. - pH) * NaturalIsotopes[3]                      +
		(float)(NEH * (NEH + 1) / 2) * powf(pH / (1. - pH), 2.) * NaturalIsotopes[2]                          -
		(float)(NEH * (NEH + 1) * (NEH + 2) / 6) * powf(pH / (1. - pH), 3.) * NaturalIsotopes[1]              +
		(float)(NEH * (NEH + 1) * (NEH + 2) * (NEH + 3) / 24) * powf(pH / (1. - pH), 4.) * NaturalIsotopes[0];

	f_rest4_XH_0 = f_rest4_XH_0 * powf(1. - pX_t / (1. - pH), NEH);

	// computing p(Rest = 4) * p(X_H = 0).

	//end of computations for I4
	//Final for I4(t) 

	fLabeledIsotopes[4] = f_rest4_XH_0 + f_rest0_XH_4 + f_rest1_XH_3 + f_rest2_XH_2 + f_rest3_XH_1;

	fLabeledIsotopes[5] = 1. - fLabeledIsotopes[0] - fLabeledIsotopes[1] - fLabeledIsotopes[2] - 
			                   fLabeledIsotopes[3] - fLabeledIsotopes[4];


	fSpectShift = Math::Abs(fLabeledIsotopes[0] - NaturalIsotopes[0]) + Math::Abs(fLabeledIsotopes[1] - NaturalIsotopes[1]);

	fSpectShift = fSpectShift + Math::Abs(fLabeledIsotopes[2] - NaturalIsotopes[2]);

	fSpectShift = fSpectShift + Math::Abs(fLabeledIsotopes[3] - NaturalIsotopes[3]);

	fSpectShift = fSpectShift + Math::Abs(fLabeledIsotopes[4] - NaturalIsotopes[4]);

	//fSpectShift = fSpectShift + Math::Abs(fLabeledIsotopes[5] - NaturalIsotopes[5]);

	//printf("SpecShift = %10.5f, fBWE = %5.4f\n", fSpectShift, fBWE1);

	//printf("Fourth ac = %10.5f \n", 
		// fLabeledIsotopes[4]);


	return 0;
}

/*
*
*   This method computes the APE of deuterium using NEH (first two mass isotopomers, fAPE1), 
*   ab initio (first three mass isotopomersm fAPE2), 
*   and statistical (first four mass isotopomers, fAPE3) 
*   
*
*/

float MassIsotopomers::ComputeAPE(array <float, 1> ^ Time_0_Isotopes, array <float, 1> ^ Time_T_Isotopes,
	float fBWE, int NEH, int Nall_Hydrogens)
{

	float fAPE1, fAPE2, fAPE2_2, fAPE3, fNEH;

	float delta_A2t_0, delta_A1t_0;

	float a, b, fRatio1, fRatio2, pXt, fdiff, fMax;

	int i, iMax;

	delta_A1t_0 = Time_T_Isotopes[1] / Time_T_Isotopes[0] - Time_0_Isotopes[1] / Time_0_Isotopes[0]; 

	delta_A2t_0 = Time_T_Isotopes[2] / Time_T_Isotopes[0] - Time_0_Isotopes[2] / Time_0_Isotopes[0];

	fAPE1 = Time_T_Isotopes[1] / Time_T_Isotopes[0] - Time_0_Isotopes[1] / Time_0_Isotopes[0];

	fAPE1 = fAPE1 * powf(1. - pH, 2.);

	fAPE1 = fAPE1 / ((float)NEH + (1. - pH) * delta_A1t_0);
		//(Time_T_Isotopes[1] / Time_T_Isotopes[0] - Time_0_Isotopes[1] / Time_0_Isotopes[0]) );
		
	

	// this part estimates enrichment, fAPE2, and number of exchnageable hydrogens, fNEH
	// from the first three mass isotopomers: M0, M1, M2.

	b = -delta_A2t_0 * (1. - pH) + Time_0_Isotopes[1] / Time_0_Isotopes[0] * pH * delta_A1t_0 * 2. * (1. - pH);

	b = b + (1. - 2. * pH) * Time_T_Isotopes[1] / Time_T_Isotopes[0] * (1. - pH) * delta_A1t_0;

	b = b + delta_A1t_0 / 2. *(delta_A1t_0 * (1. - pH) * (2. * pH - 1) + delta_A1t_0 * 2. * pH * (1. - pH) - 2. * pH);

	a = delta_A2t_0 - Time_0_Isotopes[1] / Time_0_Isotopes[0] * pH * delta_A1t_0;

	a = a - Time_T_Isotopes[1] / Time_T_Isotopes[0] * (1. - pH) * delta_A1t_0;

	a = a + delta_A1t_0 * (-delta_A1t_0 * (2. * pH - 1.) + (2.*pH - 1.) / (1. - pH)) / 2.;

	if (a != 0.0)
	{
		fAPE2 = -b / a;
	}
	
	if (fAPE2 > 0)
	{
		fNEH = delta_A1t_0 * ((1. - pH) * (1. - pH) / fAPE2 - (1. - pH));
	}
	

	// finished with fAPE1, fAPE2. Both are correct.

	//determine fAPE from the ratios of all three MIs: M0, M1, M2

	fMax = 10.0; iMax = 0;

	printf("%d\n", (int)(fBWE * 100000) );

	for (i = 0; i < (int) (fBWE* 100000); i++)
	{
		pXt = i * 0.00001;

		delta_A1t_0 = (float)NEH * pXt / ((1. - pH) * (1. - pH - pXt)) + Time_0_Isotopes[1] / Time_0_Isotopes[0];

		delta_A2t_0 = Time_0_Isotopes[2] / Time_0_Isotopes[0] - Time_0_Isotopes[1] / Time_0_Isotopes[0] * pH * (float)NEH / (1. - pH);

		delta_A2t_0 = delta_A2t_0 + powf(pH / (1. - pH), 2.) * (float)NEH * (NEH + 1) / 2;

		delta_A2t_0 = delta_A2t_0 - powf((pH + pXt)/ (1. - pH - pXt), 2.) * (float)NEH * (NEH + 1) / 2;

		delta_A2t_0 = delta_A2t_0 + (float)NEH * (pXt + pH) / (1. - pH - pXt) * delta_A1t_0;

		fdiff = 1. / (1. + delta_A1t_0 + delta_A2t_0) - Time_T_Isotopes[0] / (Time_T_Isotopes[0] + Time_T_Isotopes[1] + Time_T_Isotopes[2]);

		if (fabs(fdiff) < fMax)
		{
			fMax = fabs(fdiff);

			iMax = i;
		}

		//printf("Fdif = %10.5f %d\n", fdiff, i);
	}
	
	fAPE2_2 = (float)iMax * 0.00001;

	printf("Enrichment from TWO: %10.5f, THREE = %5.4f, and NEH = %5.1f, Three Ratio = %5.4f\n", 
		fAPE1, fAPE2, fNEH, fAPE2_2);

	return fAPE1;
}


/*
*
The most recent and accurate mass isotopomer dynamics, 09.06.2022
*
*/

int MassIsotopomers::CalculateMIDynamics(array <float, 1> ^ NaturalIsotopes,
	array <float, 1> ^ fLabeledIsotopes, float fBWE1, int NEH)
{

	float a1, a2, a3, a4, a5;

	float c1, c2, c3, c4, c5;

	float pX_t, pX_H_ratio, fSpectShift;

	float d1, d2, d3, d4;

	pX_t = fBWE1;

	pX_H_ratio = 1. - pX_t / (1. - pH);

	a1 = (float)NEH * (pX_t + pH) / (1. - pH - pX_t);

	a2 = (float) (NEH * (NEH - 1)) / 2;

	a2 = a2 * Math::Pow((pX_t + pH) / (1. - pH - pX_t), 2.0);

	a3 = (float)NEH * ((float)NEH - 1.) *  ((float)NEH - 2.) / 6.;

	a3 = a3 * Math::Pow((pX_t + pH) / (1. - pH - pX_t), 3.);

	a4 = (float) (NEH*(NEH - 1)*(NEH - 2)*(NEH-3) / 24 );   // 24 = 4!

	a4 = a4 * Math::Pow((pX_t + pH) / (1. - pH - pX_t), 4.);

	a5 = (float)(NEH*(NEH - 1)*(NEH - 2)*(NEH - 3)*(NEH - 4) / 120 );  // 120 = 5!

	a5 = a5 * Math::Pow((pX_t + pH) / (1. - pH - pX_t), 5.);

	c1 = (float)(NEH) * pH / (1. - pH);

	c2 = (float) NEH * (NEH + 1) / 2.;

	c2 = c2 * Math::Pow( pH / (1. - pH), 2.0);

	c3 = (float) NEH * (NEH + 1) * (NEH + 2) / 6;

	c3 = c3 * Math::Pow(pH / (1. - pH), 3.0);

	c4 = (float)NEH * (NEH + 1) * (NEH + 2) * (NEH + 3) / 24;

	c4 = c4 * Math::Pow(pH / (1. - pH), 4.0);

	c5 = (float)NEH * (NEH + 1) * (NEH + 2) * (NEH + 3) * (NEH + 4) / 120;

	c5 = c5 * Math::Pow(pH / (1. - pH), 5.0);



	// the monoisotope
	fLabeledIsotopes[0] = Math::Pow(pX_H_ratio, NEH) * NaturalIsotopes[0];

	fLabeledIsotopes[1] = Math::Pow(pX_H_ratio, NEH - 1) * pX_t * NEH * NaturalIsotopes[0] / Math::Pow(1. - pH, 2);

	fLabeledIsotopes[1] = fLabeledIsotopes[1] + Math::Pow(pX_H_ratio, NEH) * NaturalIsotopes[1];

	fLabeledIsotopes[2] = NaturalIsotopes[2] / NaturalIsotopes[0] - (NaturalIsotopes[1] / NaturalIsotopes[0]) * pH * (float)NEH / (1. - pH);

	fLabeledIsotopes[2] = fLabeledIsotopes[2] + Math::Pow(pH / (1. - pH), 2.) * (float)NEH * ((float)NEH + 1.) / 2.0;

	fLabeledIsotopes[2] = fLabeledIsotopes[2] - Math::Pow((pX_t + pH) / (1 - pH - pX_t), 2) * (float)NEH * ((float)NEH + 1.) / 2.0;

	fLabeledIsotopes[2] = fLabeledIsotopes[2] * fLabeledIsotopes[0];

	fLabeledIsotopes[2] = fLabeledIsotopes[2] + (float)NEH * (pX_t + pH) * fLabeledIsotopes[1] / (1. - pH - pX_t);


	//I_3(t)

	fLabeledIsotopes[3] = a3 * fLabeledIsotopes[0] + a2 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) +

		a1 * (fLabeledIsotopes[2] - a1 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) - a2 * fLabeledIsotopes[0]);


	fLabeledIsotopes[3] = fLabeledIsotopes[3] + 
		fLabeledIsotopes[0] * (NaturalIsotopes[3] / NaturalIsotopes[0] - c1 *  NaturalIsotopes[2] / NaturalIsotopes[0] +
			c2 *  NaturalIsotopes[1] / NaturalIsotopes[0] - c3);

	//I_3(t) end

	//I_4(t)

	fLabeledIsotopes[4] = a4 * fLabeledIsotopes[0] + a3 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) +

		a2 * (fLabeledIsotopes[2] - a1 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) - a2 * fLabeledIsotopes[0]);


	fLabeledIsotopes[4] = fLabeledIsotopes[4] +
		a1 * (fLabeledIsotopes[3] - a1 * (fLabeledIsotopes[2] - a1 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) - a2 * fLabeledIsotopes[0]) - 
			a2 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) - a3 * fLabeledIsotopes[0] );

	fLabeledIsotopes[4] = fLabeledIsotopes[4] +
		fLabeledIsotopes[0] * (NaturalIsotopes[4] / NaturalIsotopes[0] - c1 *  NaturalIsotopes[3] / NaturalIsotopes[0] +
			c2 *  NaturalIsotopes[2] / NaturalIsotopes[0] - c3 * NaturalIsotopes[1] / NaturalIsotopes[0] + c4);

	//I_4(t) end

	//I_5(t)

	/* fLabeledIsotopes[5] = a5 * fLabeledIsotopes[0] + (a4 - a1 * a3) * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) +

		(a3 - a1  * a2) * (fLabeledIsotopes[2] - a1 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) - a2 * fLabeledIsotopes[0]);


	 fLabeledIsotopes[5] = fLabeledIsotopes[5] +
		(a2 - a1 * a1) * (fLabeledIsotopes[3] - a1 * (fLabeledIsotopes[2] - a1 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) -
			a2 * fLabeledIsotopes[0]) - a2 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) - a3 * fLabeledIsotopes[0] );

	fLabeledIsotopes[5] = fLabeledIsotopes[5] +
		fLabeledIsotopes[0] * (NaturalIsotopes[5] / NaturalIsotopes[0] - c1 *  NaturalIsotopes[4] / NaturalIsotopes[0] +
			c2 *  NaturalIsotopes[3] / NaturalIsotopes[0] - c3 * NaturalIsotopes[2] / NaturalIsotopes[0] + 
			c4 * NaturalIsotopes[2] / NaturalIsotopes[0] - c5);
	*/

	// the firs three terms

	d1 = fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0];

	d2 = fLabeledIsotopes[2] - a1 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) - a2 * fLabeledIsotopes[0];

	fLabeledIsotopes[5] = a5 * fLabeledIsotopes[0] + a4 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) +
		a3 * (fLabeledIsotopes[2] - a1 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) - a2 * fLabeledIsotopes[0]);

	//the fourth term::

	fLabeledIsotopes[5] = fLabeledIsotopes[5] +
		a2 * (fLabeledIsotopes[3] - a1 * (fLabeledIsotopes[2] - a1 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) - a2 * fLabeledIsotopes[0]) -
			a2 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0])  - a3 * fLabeledIsotopes[0]);

	// the fifth term::

	fLabeledIsotopes[5] = fLabeledIsotopes[5] +
		a1 * (fLabeledIsotopes[4] - a4 *  fLabeledIsotopes[0]);

	d1 = fLabeledIsotopes[3] - a1 * (fLabeledIsotopes[2] - a1 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) - a2 * fLabeledIsotopes[0]) -
		a2 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) - a3 * fLabeledIsotopes[0];

	d2 = fLabeledIsotopes[2] - a1 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) - a2 * fLabeledIsotopes[0];

	d3 = fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0];

	fLabeledIsotopes[5] = fLabeledIsotopes[5] - a1 * (a1 * d1 + a2 * d2 + a3 * d3);
			

	//the sixth term
	fLabeledIsotopes[5] = fLabeledIsotopes[5] +
		fLabeledIsotopes[0] * (NaturalIsotopes[5] / NaturalIsotopes[0] - c1 *  NaturalIsotopes[4] / NaturalIsotopes[0] +
			c2 *  NaturalIsotopes[3] / NaturalIsotopes[0] - c3 * NaturalIsotopes[2] / NaturalIsotopes[0] +
			c4 * NaturalIsotopes[2] / NaturalIsotopes[0] - c5);
	printf("%10.7f\n", fLabeledIsotopes[5]);

	//I_5(t) end


	fLabeledIsotopes[5] = a5 * fLabeledIsotopes[0] + (a4 - a1 * a3) * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) +

		(a3 - a1  * a2) * (fLabeledIsotopes[2] - a1 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) - a2 * fLabeledIsotopes[0]);


	fLabeledIsotopes[5] = fLabeledIsotopes[5] + a1 * (fLabeledIsotopes[4] - a4 * fLabeledIsotopes[0]) +
		(a2 - a1 * a1) * (fLabeledIsotopes[3] - a1 * (fLabeledIsotopes[2] - a1 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) -
			a2 * fLabeledIsotopes[0]) - a2 * (fLabeledIsotopes[1] - a1 * fLabeledIsotopes[0]) - a3 * fLabeledIsotopes[0]);

	fLabeledIsotopes[5] = fLabeledIsotopes[5] +
		fLabeledIsotopes[0] * (NaturalIsotopes[5] / NaturalIsotopes[0] - c1 *  NaturalIsotopes[4] / NaturalIsotopes[0] +
			c2 *  NaturalIsotopes[3] / NaturalIsotopes[0] - c3 * NaturalIsotopes[2] / NaturalIsotopes[0] +
			c4 * NaturalIsotopes[2] / NaturalIsotopes[0] - c5);

	return 0;
}
