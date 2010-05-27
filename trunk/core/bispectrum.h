#pragma once 
#include "color.h"
#include <iostream> 
#include <fstream>
#include <sstream>

using namespace std;
class COREDLL Bispectrum : public Spectrum{
 public:
  // Spectrum Public Methods

	Bispectrum(string &filename);
	
   //perform the matrix multiplication and return the result
	Spectrum output(Spectrum & input);

 private:
	int nInputIndices;
	int nOutputIndices;
	float * inputIndices;
	float * outputIndices;
	
	float * data;

  };
