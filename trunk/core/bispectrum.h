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
	/*~Bispectrum() {
		free(inputIndices);
		free(outputIndices);
		free(data);
	}*/
	
   //perform the matrix multiplication and return the result
	Spectrum output(Spectrum & input);
	void printMyself();

 private:
	int nInputIndices;
	int nOutputIndices;
	int * inputIndices;
	int * outputIndices;
	
	float * data;

  };
