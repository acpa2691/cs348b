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
		printf("opened file '%s'\n", filename.c_str());
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
			//printf("current line of file '%s' is '%s' and has %d tokens\n", filename.c_str(), line.c_str(), (int)tokens.size());
			if(row == 0)
			{
				nOutputIndices = tokens.size();
				outputIndices = (float*)malloc(nOutputIndices*sizeof(float));
			}else if(row == 1)
			{
				nInputIndices = tokens.size();
				inputIndices = (float*)malloc(nInputIndices*sizeof(float));
				data = (float*)malloc(nInputIndices*nOutputIndices*sizeof(float));
			}
			
			for(unsigned int i = 0; i < tokens.size(); i++)
			{
				float currentFloat = atof(tokens.at(i).c_str());
				printf("current float: %f\n", currentFloat);
				if(row == 0)
				{
					outputIndices[i] = currentFloat;
				}else if(row == 1)
				{
					inputIndices[i] = currentFloat;
				}else{
					int curIndex = nOutputIndices*(row-2) + i;
					data[curIndex] = currentFloat;
				}
			}
			row++;
		}
	}else{
		printf("ERROR could not open reradiation file: %s\n", filename.c_str());
	}
	
	myfile.close();
}

Spectrum Bispectrum::output(Spectrum & input)
{
	Spectrum result(0.f); 
	
	for(int i = 0; i < nOutputIndices; i++)
	{
		//int curBaseIndex = i*nOutputIndices;
		float curOutputWavelength = outputIndices[i];
		float total = 0.f;
		for(int k = 0; k < nInputIndices; k++)
		{
			int curIndex = k*nOutputIndices+i;
			float curInputWavelength = inputIndices[k];
			total += data[curIndex]*input.getValueAtWavelength(curInputWavelength);
		}
		result.setValueAtWavelength(total, curOutputWavelength);
	}
	
	return result;
}
