// IsotopomerDynamics.h

#pragma once

using namespace System;

namespace IsotopomerDynamics {

	public ref class MassIsotopomers
	{
		// TODO: Add your methods for this class here.

	public:
		array <float> ^fINatural;
		float fBWE;

		int CalculateMIDynamics(array <float,1> ^ NaturalIsotopes, array <float, 1> ^ fLabeledIsotopes, 
			float fBWE, int NEH, int Nall_Hydrogens);

		int CalculateMIDynamics(array <float, 1> ^ NaturalIsotopes, array <float, 1> ^ fLabeledIsotopes,
			float fBWE, int NEH);

		float ComputeAPE(array <float, 1> ^ NaturalIsotopes, array <float, 1> ^ fLabeledIsotopes,
			float fBWE, int NEH, int Nall_Hydrogens);
	};
}
