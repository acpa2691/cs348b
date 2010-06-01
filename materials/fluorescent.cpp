
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

// Fluorescent.cpp*
#include "pbrt.h"
#include "material.h"
#include "reflection.h"
// Fluorescent Class Declarations
class Fluorescent : public Material {
public:
	// Fluorescent Public Methods
	Fluorescent(Reference<Texture<Spectrum> > kd,
			Reference<Texture<Spectrum> > ks,
			Reference<Texture<float> > rough,
			Reference<Texture<float> > bump,
			string & filename) {
		Kd = kd;
		Ks = ks;
		roughness = rough;
		bumpMap = bump;
		fluoroBxDF = new FluoroBxDF(filename);
	}
	BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
	              const DifferentialGeometry &dgShading) const;
private:
	// Fluorescent Private Data
	Reference<Texture<Spectrum> > Kd, Ks;
	Reference<Texture<float> > roughness, bumpMap;
	FluoroBxDF * fluoroBxDF;
};
// Fluorescent Method Definitions
BSDF *Fluorescent::GetBSDF(const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading) const {
	// Allocate _BSDF_, possibly doing bump-mapping with _bumpMap_
	DifferentialGeometry dgs;
	if (bumpMap)
		Bump(bumpMap, dgGeom, dgShading, &dgs);
	else
		dgs = dgShading;
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dgs, dgGeom.nn);
	Spectrum kd = Kd->Evaluate(dgs).Clamp();
	BxDF *diff = BSDF_ALLOC(Lambertian)(kd);
	Fresnel *fresnel =
		BSDF_ALLOC(FresnelDielectric)(1.5f, 1.f);
	Spectrum ks = Ks->Evaluate(dgs).Clamp();
	float rough = roughness->Evaluate(dgs);
	BxDF *spec = BSDF_ALLOC(Microfacet)(ks, fresnel,
		BSDF_ALLOC(Blinn)(1.f / rough));
	//bsdf->Add(diff);
	//bsdf->Add(spec);
	bsdf->Add(fluoroBxDF);
	return bsdf;
}
// Fluorescent Dynamic Creation Routine
extern "C" DLLEXPORT Material * CreateMaterial(const Transform &xform,
		const TextureParams &mp, ParamSet * paramSet = NULL) {
	printf("MAKING FLUORESCENT MATERIAL\n");
	Reference<Texture<Spectrum> > Kd = mp.GetSpectrumTexture("Kd", Spectrum(1.f));
	Reference<Texture<Spectrum> > Ks = mp.GetSpectrumTexture("Ks", Spectrum(1.f));
	Reference<Texture<float> > roughness = mp.GetFloatTexture("roughness", .1f);
	Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
	string fluoroFile = "red_ink.txt";
	if(paramSet != NULL)
	{
		fluoroFile = paramSet->FindOneString("reradiation", fluoroFile);
	}
	return new Fluorescent(Kd, Ks, roughness, bumpMap, fluoroFile);
}
