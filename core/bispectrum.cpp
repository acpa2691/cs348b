#pragma once 
#include "bispectrum.h"

Bispectrum::Bispectrum(string &filename)
{
	//parse the file
}

Spectrum Bispectrum::output(Spectrum & input)
{
	Spectrum result(0.f); 
	/*
	for(int i = 0; i < nOutputIndices; i++)
	{
		int curBaseIndex = i*nOutputIndices;
		float curOutputWavelength = outputIndices[i];
		float total = 0.f;
		for(int k = 0; k < nInputIndices; k++)
		{
			float curInputWavelength = inputIndices[k];
			total += data[curBaseIndex+k]*input.getValueAtWavelength(curInputWavelength);
		}
		result.setValueAtWavelength(total, curOutputWavelength);
	}
	*/
	return result;
}
