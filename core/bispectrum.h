#pragma once 
#include "color.h"
#include <iostream> 
#include <fstream>
#include <sstream>

using namespace std;
class COREDLL Bispectrum : public Spectrum{
 public:
  // Spectrum Public Methods

  Bispectrum(string &filename) 
  {
    //parse the file
  }
	
   //perform the matrix multiplication and return the result
   Spectrum timeSpectrum(Spectrum & input)
	{
		Spectrum result; 
		
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
		
		return result;
	}

 private:
	int nInputIndices;
	int nOutputIndices
	float * inputIndices;
	float * outputIndices;
	
	float * data;

  };
