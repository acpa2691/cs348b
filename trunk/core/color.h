

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
#include <algorithm>

#define SPECTRUM_SAMPLES nCIE
#define SPECTRUM_SPACING 1
#define SPECTRUM_START CIEstart
#define SPECTRUM_END CIEend

using namespace std;

// Spectrum Declarations
class COREDLL Spectrum {
 public:
  // Spectrum Public Methods

  Spectrum(float v = 0.f) {
    fill(c,c+SPECTRUM_SAMPLES,v);
  }
	//cs needs to have sive SPECTRUM_SAMPLES
  Spectrum(float * cs) {
    for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
      c[i] = cs[i];
  }
  
Spectrum(float mean, float stdev, float scale){
    fill(c,c+SPECTRUM_SAMPLES,0.f);
    int minLambda = max(mean - 3 * stdev, (float)SPECTRUM_START);
    int maxLambda = min (mean + 3 * stdev, (float)SPECTRUM_END);
    int index = minLambda - SPECTRUM_START;
    float var = stdev*stdev;
    float coeff = scale / sqrt(2*3.145*var);
    float invTwoVar = -1.f / (2.f*var);
	printf("index:%d minLambda: %d maxLambda: %d\n", index, minLambda, maxLambda);
    while(index <= maxLambda - SPECTRUM_START){
      float lambda = index + SPECTRUM_START;
      c[index] = coeff * expf(invTwoVar * (lambda - mean)*(lambda-mean));
		printf("added value: %f\n", c[index]);
		index++;
    }
  }

    friend ostream &operator<<(ostream &, const Spectrum &);
    Spectrum &operator+=(const Spectrum &s2) {
		//printf("attempting add\n");
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	{
		c[i] += s2.c[i];
		//printf("adding %d\n", i);
      }
	return *this;
    }
    Spectrum operator+(const Spectrum &s2) const {
		//printf("attempting add\n");
      Spectrum ret = *this;
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	  {
		  ret.c[i] += s2.c[i];
		 // printf("adding %d\n", i);
	  }
      return ret;
    }
    Spectrum operator-(const Spectrum &s2) const {
      Spectrum ret = *this;
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	ret.c[i] -= s2.c[i];
      return ret;
    }
    Spectrum operator/(const Spectrum &s2) const {
      Spectrum ret = *this;
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	ret.c[i] /= s2.c[i];
      return ret;
    }
    Spectrum operator*(const Spectrum &sp) const {
      Spectrum ret = *this;
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	ret.c[i] *= sp.c[i];
      return ret;
    }
    Spectrum &operator*=(const Spectrum &sp) {
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	c[i] *= sp.c[i];
      return *this;
    }
    Spectrum operator*(float a) const {
      Spectrum ret = *this;
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	ret.c[i] *= a;
      return ret;
    }
    Spectrum &operator*=(float a) {
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	c[i] *= a;
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
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	c[i] *= inv;
      return *this;
    }
    void AddWeighted(float w, const Spectrum &s) {
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	c[i] += w * s.c[i];
    }
    bool operator==(const Spectrum &sp) const {
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	if (c[i] != sp.c[i]) return false;
      return true;
    }
    bool operator!=(const Spectrum &sp) const {
      return !(*this == sp);
    }
    bool Black() const {
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	if (c[i] != 0.) return false;
      return true;
    }
    Spectrum Sqrt() const {
      Spectrum ret;
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	ret.c[i] = sqrtf(c[i]);
      return ret;
    }
    Spectrum Pow(const Spectrum &e) const {
      Spectrum ret;
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	ret.c[i] = c[i] > 0 ? powf(c[i], e.c[i]) : 0.f;
      return ret;
    }
    Spectrum operator-() const {
      Spectrum ret;
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	ret.c[i] = -c[i];
      return ret;
    }
    friend Spectrum Exp(const Spectrum &s) {
      Spectrum ret;
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	ret.c[i] = expf(s.c[i]);
      return ret;
    }
    Spectrum Clamp(float low = 0.f,
		   float high = INFINITY) const {
      Spectrum ret;
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	ret.c[i] = ::Clamp(c[i], low, high);
      return ret;
    }
    bool IsNaN() const {
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	if (isnan(c[i])) return true;
      return false;
    }
    void Print(FILE *f) const {
      for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
	fprintf(f, "%f ", c[i]);
    }
    void XYZ(float xyz[3]) const {
      xyz[0] = xyz[1] = xyz[2] = 0.;
      for (int i = SPECTRUM_START; i < SPECTRUM_END; ++i) {
	if( i > CIEstart && i < CIEend){
	  xyz[0] += XWeight[i-CIEstart] * c[i-SPECTRUM_START];
	  xyz[1] += YWeight[i-CIEstart] * c[i-SPECTRUM_START];
	  xyz[2] += ZWeight[i-CIEstart] * c[i-SPECTRUM_START];
	}
      }
    }
    float y() const {
      float v = 0.;
      for (int i = SPECTRUM_START; i < SPECTRUM_END; ++i) 
	  {
		if( i > CIEstart && i < CIEend){
		  v += YWeight[i-CIEstart] * c[i-SPECTRUM_START];
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
    float c[SPECTRUM_SAMPLES];
    static float XWeight[COLOR_SAMPLES];
    static float YWeight[COLOR_SAMPLES];
    static float ZWeight[COLOR_SAMPLES];
    friend Spectrum FromXYZ(float x, float y, float z);
  };
#endif // PBRT_COLOR_H
