

/*
 pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys.
 
 This file is part of pbrt.
 
 pbrt is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.  Note that the text contents of
 the book "Physically Based Rendering" are *not* licensed under the
 GNU GPL.
 
 pbrt is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 */

#ifndef PBRT_COLOR_H
#define PBRT_COLOR_H
// color.h*
#include "pbrt.h"
#include <map>
#include <algorithm>

//#define SPECTRUM_SAMPLES 471
//#define SPECTRUM_SPACING 1
//#define SPECTRUM_START CIEstart
//#define SPECTRUM_END SPECTRUM_START + SPECTRUM_SAMPLES
//#define SPECTRUM_SAMPLES (SPECTRUM_END - SPECTRUM_START + 1)

#define WAVELENGTH_RED 575
#define WAVELENGTH_GREEN 535
#define WAVELENGTH_BLUE 445

#define VISIBLE_START 400
#define VISIBLE_END 700

using namespace std;
typedef map<int,float>::iterator sample_iterator;
typedef map<int,float>::const_iterator const_sample_iterator;

// Spectrum Declarations
class COREDLL Spectrum {
public:
	// Spectrum Public Methods
	
	Spectrum(float v = 0.f) {
		onlyVisible = true;
		defaultScale = v;
		/*
		 for(int i = 0; i < nCIE; i++)
		 {
		 samples[CIEstart+i] = defaultScale;
		 }
		 */
		samples[WAVELENGTH_BLUE]= defaultScale;
		samples[WAVELENGTH_GREEN] = defaultScale;
		samples[WAVELENGTH_RED] = defaultScale;
	}
	//cs needs to have sive SPECTRUM_SAMPLES
	Spectrum(float  cs[3]) {
		onlyVisible = true;
		defaultScale = 0.f;
		samples[WAVELENGTH_BLUE]= cs[2];
		samples[WAVELENGTH_GREEN] = cs[1];
		samples[WAVELENGTH_RED] = cs[0];
		if(cs[0] == cs[1] && cs[1] == cs[2])
		{
			defaultScale = cs[0];
		}
		/*
		 for(int i = 0; i < nCIE; i++)
		 {
		 samples[CIEstart+i] = cs[0];
		 }
		 */
	}
	
	void addGaussian(float mean, float stdev, float height)
	{
		int minLambda = (int)(mean - 3 * stdev);
		int maxLambda = (int)(mean + 3* stdev);
		float var = stdev*stdev;
		float coeff = height / sqrt(2*3.145*var);
		float invTwoVar = -1.f / (2.f*var);
		float step = 1.f;
		float spreadFactor = 1.3f;
		int nSamples = 1;
		
		setValueAtWavelength(coeff, mean);
		while(true)
		{
			if(mean - step < minLambda || mean + step > maxLambda)
			{
				break;
			}
			int curStep = (int)step;
			float curValue = coeff*expf(invTwoVar*curStep*curStep);
			setValueAtWavelength(curValue, (int)(mean + curStep));
			setValueAtWavelength(curValue, (int)(mean - curStep));
			step += spreadFactor*coeff/curValue;
			nSamples += 2;
			//printf("adding value %f at step %d\n", curValue, mean+curStep);
		}
		
		//printf("made gaussian with %d samples ", nSamples);
		//printSelf();
		
		/*
		 //printf("creating spectrum with var: %f coeff: %f invTwoVar: %f\n", var, coeff, invTwoVar);
		 for(int lambda = minLambda; lambda < maxLambda; lambda ++){
		 float curValue = coeff * expf(invTwoVar * (lambda - mean)*(lambda-mean));
		 //printf("adding value: %f for lambda: %d\n", curValue, lambda);
		 setValueAtWavelength(curValue, lambda);
		 }
		 */
		
	}
	
	Spectrum(float mean, float stdev, float scale){
		defaultScale = 0.f;
		addGaussian(mean, stdev, scale);
		//int minLambda = max(mean - 3 * stdev, (float)SPECTRUM_START);
		//int maxLambda = min (mean + 3 * stdev, (float)SPECTRUM_END);
		
		//printf("and now printing myself: ");
		//printSelf();
	}
	
	static bool SpectrumTest();
	void printSelf() const;
	
	friend ostream &operator<<(ostream &, const Spectrum &);
	
	void setValueAtWavelength(float val, int lambda)
	{
		if((lambda < VISIBLE_START || lambda > VISIBLE_END) && val > 0)
		{
			onlyVisible = false;
		}
		samples[lambda] = val;
	}
	
	float getValueAtWavelength(int lambda) const
	{
		if(onlyVisible)
		{
			if(lambda < VISIBLE_START || lambda > VISIBLE_END)
			{
				return 0.f;
			}
		}
		const_sample_iterator itr  = samples.find(lambda);
		if(itr != samples.end()){
			return itr->second;
		}
		else{
			const_sample_iterator before = samples.lower_bound(lambda);
			before--;
			const_sample_iterator after = samples.upper_bound(lambda);
			
			if(before == samples.begin() || after == samples.end() || before == after)
				return defaultScale;
			
			if(before == after)
			{
				return before->second;
			}
			
			return (lambda - before->first) * (after->second - before->second)/(after->first - before->first);
		}
	}
	
	Spectrum &operator+=(const Spectrum &s2) {
		for (const_sample_iterator itr = s2.samples.begin(); itr != s2.samples.end();
			 ++itr){
			setValueAtWavelength(itr->second + getValueAtWavelength(itr->first), itr->first);
		}
		return *this;
	}
	Spectrum operator+(const Spectrum &s2) const {
		//printf("attempting add\n");
		Spectrum ret = *this;
		for (const_sample_iterator itr = s2.samples.begin(); itr != s2.samples.end();
			 ++itr){
			ret.setValueAtWavelength(itr->second + ret.getValueAtWavelength(itr->first), itr->first);
		}
		return ret;
	}
	Spectrum operator-(const Spectrum &s2) const {
		Spectrum ret = *this;
		for (const_sample_iterator itr = s2.samples.begin(); itr != s2.samples.end();
			 ++itr){
			ret.setValueAtWavelength(ret.getValueAtWavelength(itr->first) - itr->second, itr->first);
		}
		return ret;
	}
	Spectrum operator/(const Spectrum &s2) const {
		Spectrum ret = *this;
		for (const_sample_iterator itr = s2.samples.begin(); itr != s2.samples.end();
			 ++itr){
			ret.setValueAtWavelength(ret.getValueAtWavelength(itr->first)/itr->second, itr->first);
		}
		return ret;
	}
	Spectrum operator*(const Spectrum &sp) const {
		Spectrum ret = *this;
		for (const_sample_iterator itr = sp.samples.begin(); itr != sp.samples.end();
			 ++itr){
			ret.setValueAtWavelength(ret.getValueAtWavelength(itr->first)*itr->second, itr->first);
		}
		return ret;
	}
	Spectrum &operator*=(const Spectrum &sp) {
		for (const_sample_iterator itr = sp.samples.begin(); itr != sp.samples.end();
			 ++itr){
			setValueAtWavelength(getValueAtWavelength(itr->first)*itr->second, itr->first);
		}
		return *this;
	}
	Spectrum operator*(float a) const {
		Spectrum ret = *this;
		for (const_sample_iterator itr = samples.begin(); itr != samples.end();
			 ++itr){
			ret.setValueAtWavelength(this->getValueAtWavelength(itr->first)*a, itr->first);
		}
		return ret;
	}
	Spectrum &operator*=(float a) {
		for (sample_iterator itr = samples.begin(); itr != samples.end();
			 ++itr){
			setValueAtWavelength(getValueAtWavelength(itr->first)*a, itr->first);
		}
		return *this;
	}
	friend inline
    Spectrum operator*(float a, const Spectrum &s) {
		return s * a;
	}
	Spectrum operator/(float a) const {
		return *this * (1.f / a);
	}
	Spectrum &operator/=(float a) {
		float inv = 1.f / a;
		for (sample_iterator itr = samples.begin(); itr != samples.end();
			 ++itr){
			setValueAtWavelength(getValueAtWavelength(itr->first)* inv, itr->first);
		}
		return *this;
	}
	void AddWeighted(float w, const Spectrum &s) {
		//printf("before addWeighted: ");
		//printSelf();
		for (const_sample_iterator itr = s.samples.begin(); itr != s.samples.end(); ++itr){
			setValueAtWavelength(s.getValueAtWavelength(itr->first) * w + getValueAtWavelength(itr->first), itr->first);
		}
		
		//printf("after addWeighted: ");
		//printSelf();
	}
	bool operator==(const Spectrum &sp) const {
		const_sample_iterator itr1, itr2;
		for(itr1 = samples.begin(), itr2 = sp.samples.begin();
			itr1 != samples.end() && itr2 != sp.samples.end(); ++itr1,
			++itr2){
			if(itr1->first != itr2 ->first || itr1->second != itr2->second)
				return false;
		}
		if(itr1 == samples.end() && itr2 == sp.samples.end()) return true;
		return false;
	}
	bool operator!=(const Spectrum &sp) const {
		return !(*this == sp);
	}
	bool Black() const {
		for (const_sample_iterator itr = samples.begin(); itr != samples.end(); ++itr)
			if (itr->second != 0. || defaultScale != 0.f) 
				return false;
		return true;
	}
	Spectrum Sqrt() const {
		Spectrum ret;
		for (const_sample_iterator itr = samples.begin(); itr !=
			 samples.end(); ++itr){
			int lambda = itr->first;
			ret.setValueAtWavelength(sqrtf(samples.find(lambda)->second), lambda);
		}
		return ret;
	}
	Spectrum Pow(const Spectrum &e) const {
		Spectrum ret;
		for (const_sample_iterator itr = e.samples.begin(); itr != e.samples.end();
			 ++itr){
			float lambda = itr->first;
			float curVal = getValueAtWavelength(lambda);
			ret.setValueAtWavelength(curVal > 0.f ? curVal : 0.f, lambda);
		}
		return ret;
	}
	Spectrum operator-() const {
		Spectrum ret;
		for (const_sample_iterator itr = samples.begin(); itr != samples.end();
			 ++itr){
			ret.setValueAtWavelength(-ret.getValueAtWavelength(itr->first), itr->first);
		}
		return ret;
	}
	
	friend Spectrum Exp(const Spectrum &s) {
		Spectrum ret;
		for (const_sample_iterator itr = s.samples.begin(); itr != s.samples.end(); ++itr)
			ret.setValueAtWavelength(expf(itr->second), itr->first);
		return ret;
	}
	Spectrum Clamp(float low = 0.f,
				   float high = INFINITY) const {
		Spectrum ret;
		for (const_sample_iterator itr = samples.begin(); itr != samples.end(); ++itr)
			ret.setValueAtWavelength(::Clamp(itr->second,low,high), itr->first);
		return ret;
	}
	bool IsNaN() const {
		for (const_sample_iterator itr = samples.begin(); itr != samples.end();   ++itr)
			if(isnan(itr->second)) return true;
		return false;
	}
	void Print(FILE *f) const {
		for (const_sample_iterator itr = samples.begin(); itr != samples.end();   ++itr)
			fprintf(f, "%f ", itr->second);
	}
	
	void XYZ(float xyz[3]) const {
		xyz[0] = xyz[1] = xyz[2] = 0.;
		//printf("XYZ printing:");
		for (const_sample_iterator itr = samples.begin(); itr != samples.end();   ++itr){
			if( itr->first > CIEstart && itr->first < CIEend){
				//printf(" (%d, %f)", itr->first, itr->second);
				float curVal = itr->second;
				int curLambda = itr->first;
				//printf(" cie x:%f y:%f z:%f ", CIE_X[curLambda - CIEstart],CIE_Y[curLambda - CIEstart],CIE_Z[curLambda - CIEstart]);
				xyz[0] += CIE_X[curLambda - CIEstart] * curVal ;
				xyz[1] += CIE_Y[curLambda - CIEstart] * curVal ;
				xyz[2] += CIE_Z[curLambda - CIEstart] * curVal ;
			}
		}
		// printf("\n");
		//printf("XYZ x:%f y:%f z:%f \n", xyz[0], xyz[1], xyz[2]);
	}
	
	float y() const {
		//printf("getting Y value for spectrum: ");
		//printSelf();
		float v = 0.;
		for (const_sample_iterator itr = samples.begin(); itr != samples.end();   ++itr){
			if( itr->first > CIEstart && itr->first < CIEend){
				float curVal = defaultScale;
				if(!samples.empty())
					curVal = itr->second;
				v += CIE_Y[itr->first - CIEstart] * curVal;
			}
		}
		return v;
	}
	bool operator<(const Spectrum &s2) const {
		return y() < s2.y();
	}
	friend class ParamSet;
	
	// Spectrum Public Data
	static const int CIEstart = 360;
	static const int CIEend = 830;
	static const int nCIE = CIEend-CIEstart+1;
	static const float CIE_X[nCIE];
	static const float CIE_Y[nCIE];
	static const float CIE_Z[nCIE];
private:
	// Spectrum Private Data
	//float c[SPECTRUM_SAMPLES];
	map<int, float> samples;
	float defaultScale;
	bool onlyVisible;
	static float XWeight[COLOR_SAMPLES];
	static float YWeight[COLOR_SAMPLES];
	static float ZWeight[COLOR_SAMPLES];
	friend Spectrum FromXYZ(float x, float y, float z);
};
#endif // PBRT_COLOR_H
