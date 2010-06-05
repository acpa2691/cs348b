	#include "bispectrum.h"
#include <iostream>
#include <fstream>
#include <string>

void Tokenize(const string& str,
			  vector<string>& tokens,
			  const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);
	
    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

Bispectrum::Bispectrum(string &filename)
{
	
	ifstream myfile;
	myfile.open(filename.c_str());
	if (myfile.is_open())
	{
		//printf("opened file '%s'\n", filename.c_str());
		int row = 0;
		string line;
		while (!myfile.eof())
		{
			getline(myfile, line);
			vector<string> tokens;
			
			Tokenize(line, tokens);
			if(tokens.size() == 0)
			{
				break;
			}
			//printf("current line(%d) of file '%s' is '%s' and has %d tokens\n", row, filename.c_str(), line.c_str(), (int)tokens.size());
			if(row == 0)
			{
				nOutputIndices = (int)tokens.size();
				outputIndices = (int*)malloc(nOutputIndices*sizeof(int));
			}else if(row == 1)
			{
				nInputIndices = (int)tokens.size();
				inputIndices = (int*)malloc(nInputIndices*sizeof(int));
				data = (float*)malloc(nInputIndices*nOutputIndices*sizeof(float));
			}
			
			for(unsigned int i = 0; i < tokens.size(); i++)
			{
				float currentFloat = atof(tokens.at(i).c_str());
				//printf("current float: %f\n", currentFloat);
				if(row == 0)
				{
					outputIndices[i] = (int)currentFloat;
				}else if(row == 1)
				{
					inputIndices[i] = (int)currentFloat;
				}else{
					int curIndex = nOutputIndices*(row-2) + i;
					data[curIndex] = 0.01*currentFloat;
				}
			}
			row++;
		}
	}else{
		printf("ERROR could not open reradiation file: %s\n", filename.c_str());
	}
	
	myfile.close();
}

void Bispectrum::scale(float scaleFactor)
{
	int total = nOutputIndices*nInputIndices;
	
	for(int i = 0; i < total; i++)
	{
		data[i] = data[i]*scaleFactor;
	}
}

void Bispectrum::printMyself()
{
	printf("i have nInputs: %d nOutputs: %d\n", nInputIndices, nOutputIndices);
}

Spectrum Bispectrum::output(Spectrum & input)
{
	return Bispectrum::output(input, true, true);
}

Spectrum Bispectrum::output(Spectrum & input, bool mainDiag, bool reemission)
{
	//printf("INPUT ");
	//input.printSelf();
	Spectrum result(0.f); 
	//printf("nInputs: %d nOutputs:%d\n", nInputIndices, nOutputIndices);
	
	for(int i = 0; i < nInputIndices; i++)
	{
		int curBaseIndex = i*nOutputIndices;
		int curInputWavelength = inputIndices[i];
		float curInputValue = input.getValueAtWavelength(curInputWavelength);
		//printf("at input #%d with wavelength: %d and value: %f\n", i+1, curInputWavelength, curInputValue);
		Spectrum currentSpec(0.f);
		for(int k = 0; k < nOutputIndices; k++)
		{
			int curIndex = curBaseIndex+k;
			int curOutputWavelength = outputIndices[k];
			//printf("at index: %d\n", curIndex);
			
			if(i == k)
			{
				if(mainDiag)
				{
					currentSpec.setValueAtWavelength(data[curIndex] * curInputValue, curOutputWavelength);
				}else{
					currentSpec.setValueAtWavelength(0.f, curOutputWavelength);
				}
			}else{
				if(reemission)
				{
					currentSpec.setValueAtWavelength(data[curIndex] * curInputValue, curOutputWavelength);
				}else{
					currentSpec.setValueAtWavelength(0.f, curOutputWavelength);
				}
			}
		}
		result += currentSpec;
	}
	
	//printf("OUTPUTTING ");
	//result.printSelf();
	
	return result;
}
