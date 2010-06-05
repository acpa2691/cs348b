
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

// path.cpp*
#include "pbrt.h"
#include "transport.h"
#include "scene.h"
#include "color.h"
#include "mc.h"
#include "primitive.h"
#include "octree.h"
#include "material.h"
#include "../materials/fluorescent.cpp"
#include "../shapes/trianglemesh.cpp"
#include <iostream>
using namespace std;


// DirectLighting Declarations
enum LightStrategy { SAMPLE_ALL_UNIFORM, SAMPLE_ONE_UNIFORM,
	SAMPLE_ONE_WEIGHTED  // NOBOOK
};

// PathIntegrator Declarations
class PathSubsurface : public SurfaceIntegrator {
public:
  // PathSubsurface Public Methods
  Spectrum Li(const Scene *scene, const RayDifferential &ray, const Sample *sample, float *alpha) const;
  void RequestSamples(Sample *sample, const Scene *scene);
  void Preprocess(const Scene * scene);
  PathSubsurface(int md) { maxDepth = md; }
private:
  // PathSubsurface Private Data
  int maxDepth;
#define SAMPLE_DEPTH 3
  int lightPositionOffset[SAMPLE_DEPTH];
  int lightNumOffset[SAMPLE_DEPTH];
  int bsdfDirectionOffset[SAMPLE_DEPTH];
  int bsdfComponentOffset[SAMPLE_DEPTH];
  int outgoingDirectionOffset[SAMPLE_DEPTH];
  int outgoingComponentOffset[SAMPLE_DEPTH];

  // DirectLighting Private Data
  LightStrategy strategy;
  mutable int rayDepth; // NOBOOK
  // Declare sample parameters for light source sampling
  //int *lightSampleOffset, lightNumOffset;
  //mutable float *avgY, *avgYsample, *cdf;
  //mutable float overallAvgY;
  
  
  float scatteringCoefficient;				// (sigma_s)
  float absorptionCoefficient; 				// (sigma_a)
  float extinctionCoefficient; 				// (sigma_t)
  float reducedExtinctionCoefficient;			// (sigma_t')
  float effectiveExtinctionCoefficient;		// (sigma_tr)
  float relativeRefractiveIndex;				// (eta)
  float albedo;								// (alpha)
  float diffusionConstant; 					// (D)
  float diffuseReflectionCoefficent;			// fraction of light reflected
  float fresnelReflectance;					// (Fr)
  float A;  									// Fresnel factor
  float g;									// mean cosine of scattering angle
};
// PathSubsurface Method Definitions
void PathSubsurface::RequestSamples(Sample *sample,
				    const Scene *scene) {
  for (int i = 0; i < SAMPLE_DEPTH; ++i) {
    lightPositionOffset[i] = sample->Add2D(1);
    lightNumOffset[i] = sample->Add1D(1);
    bsdfDirectionOffset[i] = sample->Add2D(1);
    bsdfComponentOffset[i] = sample->Add1D(1);
    outgoingDirectionOffset[i] = sample->Add2D(1);
    outgoingComponentOffset[i] = sample->Add1D(1);
  }
}
Spectrum PathSubsurface::Li(const Scene *scene,
			    const RayDifferential &r, const Sample *sample,
			    float *alpha) const {
  // Declare common path integration variables
  Spectrum pathThroughput = 1., L = 0.;
  RayDifferential ray(r);
  bool specularBounce = false;

  for (int pathLength = 0; ; ++pathLength) {
    cout <<"path length " << pathLength <<endl;
    // Find next vertex of path
    Intersection isect;
    if (!scene->Intersect(ray, &isect)) {
      // Stop path sampling since no intersection was found
      if (pathLength == 0) {
	for (u_int i = 0; i < scene->lights.size(); ++i)
	  L += scene->lights[i]->Le(ray);
	if (alpha) {
	  if (L != 0.) *alpha = 1.;
	  else *alpha = 0.;
	}
      }
      else if (pathLength > 0 && alpha && specularBounce) {
	for (u_int i = 0; i < scene->lights.size(); ++i)
	  L += pathThroughput * scene->lights[i]->Le(ray);
      }
      break;
    }
    if (pathLength == 0) {
      r.maxt = ray.maxt;
      if (alpha) *alpha = 1.;
    }
    else
      {
	pathThroughput *= scene->Transmittance(ray);
      }
    // Possibly add emitted light at path vertex
    if (pathLength == 0 || specularBounce)
      {
	L += pathThroughput * isect.Le(-ray.d);
      }

    BSDF *bsdf = isect.GetBSDF(ray);
    // Sample illumination from lights to find path contribution
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    Vector wo = -ray.d;
		
		
    if (scene->lights.size() > 0) {
      // Apply subsurface scattering if the material is fluorescent
      if (isect.primitive->IsFluorescent()) {
	//printf("at fluorescent intersection\n");
	const GeometricPrimitive* prim = (GeometricPrimitive*)isect.primitive;
	Reference<Shape> shapeRef = prim->shape;
	TriangleMesh * theShape = (TriangleMesh*)(shapeRef.getObject());
	int cacheIndex = 0;//theShape->CacheIndex();
	vector<irradianceCache>& cache = FluorescentShapesCache[cacheIndex];
	//printf("at cache index:%d with size: %d\n", cacheIndex, (int)cache.size());
	assert(prim->material->IsFluorescent());
				
				
	Point& xo = isect.dg.p;
	//printf("abotu to get material\n");
	const FluorescentMaterial * material = (FluorescentMaterial*)(prim->material.getObject());
	//printf("got material\n");
	float lu = material->meanFreePath;
	//printf("accessed material\n");
	float zv = material->virtualLightSourceHeight;
	float A = material->A;
	Bispectrum * fluoro = (Bispectrum*)(bsdf->f_ptr(wo, wo, BxDFType(BSDF_FLUORESCENT)));
	Spectrum sigmatr = material->effectiveTransportCoefficient;
	Spectrum alphaprime = material->reducedAlbedo;
	//printf("FULLY accessed material\n");
	for (int i = 0; i < (int)cache.size(); i++) {
	  //printf("going through cache\n");
	  Spectrum& E = cache[i].E;
	  if(!E.Black())
	    {
	      Point& xi = cache[i].p;
	      Vector n(cache[i].n);
	      float Area = cache[i].Area;
	      //printf("i have cached irradiance ");
	      //E.printSelf();
	      Point Pr = xi - lu * n;
	      float dr = (xo - Pr).Length();
	      L +=  (1.f/(dr*dr))* Area * E; // * Ft
	    }
	}
	//printf("this is radiance BEFORE ");
	//L.printSelf();
	L *= AbsDot(wo, n);
	L = pathThroughput * fluoro->output(L);
	//printf("this is radiance AFTER ");
	//L.printSelf();
      }
      else{
	if (pathLength < SAMPLE_DEPTH)
	  {
	    L += pathThroughput * UniformSampleOneLight(scene, p, n,wo, bsdf, sample,lightPositionOffset[pathLength],lightNumOffset[pathLength], bsdfDirectionOffset[pathLength], bsdfComponentOffset[pathLength]);
	  }
	else
	  {
	    L += pathThroughput * UniformSampleOneLight(scene, p, n, wo, bsdf, sample);
	  }
      }
		
		
      // Sample BSDF to get new path direction
      // Get random numbers for sampling new direction, _bs1_, _bs2_, and _bcs_
      float bs1, bs2, bcs;
      if (pathLength < SAMPLE_DEPTH) {
	bs1 = sample->twoD[outgoingDirectionOffset[pathLength]][0];
	bs2 = sample->twoD[outgoingDirectionOffset[pathLength]][1];
	bcs = sample->oneD[outgoingComponentOffset[pathLength]][0];
      }
      else {
	bs1 = RandomFloat();
	bs2 = RandomFloat();
	bcs = RandomFloat();
      }
      Vector wi;
      float pdf;
      BxDFType flags;


      if(bsdf->NumComponents(BSDF_FLUORESCENT) > 0)
	{
	  Bispectrum * fluoro = (Bispectrum*)bsdf->Sample_f_ptr(wo, &wi, bs1, bs2, bcs, &pdf, BSDF_FLUORESCENT, &flags);
	  pathThroughput = fluoro->output(pathThroughput)* AbsDot(wi, n) / pdf;
	  specularBounce = false;
	  if(pdf == 0)
	    {
	      cout << "PDF is O "<<endl;
	      break;
	    }
	}	
	       
		
      Spectrum f = bsdf->Sample_f(wo, &wi, bs1, bs2, bcs, &pdf, BxDFType(BSDF_ALL & ~BSDF_FLUORESCENT), &flags);
      if (f.Black() || pdf == 0.){
	cout <<"f.Black() || pdf == 0." <<endl;
	break;
      }
      specularBounce = (flags & BSDF_SPECULAR) != 0;
      pathThroughput *=  f * AbsDot(wi, n) / pdf;
		
      ray = RayDifferential(p, wi);
      //Possibly terminate the path
      if (pathLength > 3) {
	float continueProbability = .5f;
	if (RandomFloat() > continueProbability)
	  break;
	pathThroughput /= continueProbability;
      }
		
       if (pathLength == maxDepth){
	 cout <<"pathLength == maxDepth"<<endl;
	break;
       }

    }
  }
  return L;
}

void PathSubsurface::Preprocess(const Scene * scene)
{
  printf("SUBSURFACE preprocessing\n");
  for (int i = 0; i < (int)FluorescentShapes.size(); i++) {
    vector<Reference<Shape> > refined;
    FluorescentShapes[i]->Refine(refined); // will generate samples
  }
  for (int i = 0; i < (int)FluorescentShapesCache.size(); i++) {
    vector<irradianceCache>& cache = FluorescentShapesCache[i];
    for (int j = 0; j < (int)cache.size(); j++) {
      irradianceCache& cacheEntry = cache[j];
      cacheEntry.E = EstimateIrradiance(scene, cacheEntry.p, cacheEntry.n, NULL);
    }
  }
}

 extern "C" DLLEXPORT SurfaceIntegrator *CreateSurfaceIntegrator(const
								 ParamSet &params) 
 {
   int maxDepth = params.FindOneInt("maxdepth", 5);
   return new PathSubsurface(maxDepth);
 }
