
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

//Modified from [song,yu]
// FluorescentMaterial Class Declarations
class  FluorescentMaterial : public Material {
public:
	FluorescentMaterial(Reference<Texture<Spectrum> > kd,
					Reference<Texture<Spectrum> > ks,
					Reference<Texture<float> > rough,
					Reference<Texture<float> > bump,
					string & filename,
					const Spectrum& reducedScatteringCoefficient,
					const Spectrum& absorptionCoefficient,
					float relativeRefractiveIndex);
	void Initialize(const Spectrum& reducedScatteringCoefficient,
		const Spectrum& absorptionCoefficient,
		float relativeRefractiveIndex);

	// FluorescentMaterial Interface
	virtual BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const;
	virtual ~FluorescentMaterial();
	virtual bool IsFluorescent() const { return true; }

	Spectrum scatteringCoefficient;				// (sigma_s)
	Spectrum reducedScatteringCoefficient;			// (sigma_s')
	Spectrum absorptionCoefficient;				// (sigma_a)
	Spectrum extinctionCoefficient;				// (sigma_t)
	Spectrum reducedExtinctionCoefficient;			// (sigma_t')
	Spectrum effectiveTransportCoefficient;			// (sigma_tr)
	//Spectrum diffuseReflectionCoefficient;		// fraction of light reflected
		
	float relativeRefractiveIndex;			// (eta)
	float diffuseFresnelReflectance;			// (Fdr)
	float diffuseFresnelTransmittance;		// (Ftr = 1 - Fdr)
	float A;  							// Fresnel factor
		
	Spectrum reducedAlbedo;
	
	float albedo;								// (alpha)
	float diffusionConstant; 					// (D)
	float fresnelReflectance;					// (Fr)
	float meanFreePath;
	float virtualLightSourceHeight;				// (zv)
	
	float g;									// mean cosine of scattering angle	
private:
	// Fluorescent Private Data
	Reference<Texture<Spectrum> > Kd, Ks;
	Reference<Texture<float> > roughness, bumpMap;
	FluoroBxDF * fluoroBxDF;
};

//Modified from [song,yu]
/************************
 * Fluorescent Material *
 ************************/
 

FluorescentMaterial::FluorescentMaterial(Reference<Texture<Spectrum> > kd, Reference<Texture<Spectrum> > ks, Reference<Texture<float> > rough, Reference<Texture<float> > bump, string & filename, const Spectrum& reducedScatteringCoefficient,
		const Spectrum& absorptionCoefficient,
		float relativeRefractiveIndex)
{
	Kd = kd;
	Ks = ks;
	roughness = rough;
	bumpMap = bump;
	fluoroBxDF = new FluoroBxDF(filename);
	
	Initialize(reducedScatteringCoefficient, absorptionCoefficient, relativeRefractiveIndex);
}


void FluorescentMaterial::Initialize(const Spectrum& reducedScatteringCoefficient,
		const Spectrum& absorptionCoefficient,
		float relativeRefractiveIndex)
{
	this->reducedScatteringCoefficient = reducedScatteringCoefficient;
	this->absorptionCoefficient = absorptionCoefficient;
	this->relativeRefractiveIndex = relativeRefractiveIndex;

	float eta = relativeRefractiveIndex;
	float invEta = 1.0/eta;
	// FresnelDielectric Fr(1.0f, eta);
	// Fr.Evaluate(cosi);
	// Fr.Evaluate(coso);
	this->diffuseFresnelReflectance = -1.440*invEta*invEta + 0.710*invEta + 0.668 + 0.0636*eta;
	this->diffuseFresnelTransmittance = 1.0 - diffuseFresnelReflectance; // by conservation of energy since dialectrics don't absorb light
	float Fdr = diffuseFresnelReflectance;
	this->A = (1.0 + Fdr) / (1.0 - Fdr);
	
	this->reducedExtinctionCoefficient = reducedScatteringCoefficient + absorptionCoefficient;  // sigma_t' = sigma_s' + sigma_a
	this->meanFreePath = 1.0 / reducedExtinctionCoefficient.y();								// l_u      = 1 / sigma_t'
	this->virtualLightSourceHeight = meanFreePath * (1 + (4.0/3.0)*A);							// zv = lu * ( 1 + 4A/3 )
	this->reducedAlbedo = reducedScatteringCoefficient / reducedExtinctionCoefficient;          // alpha'   = sigma_s' / sigma_t'
	this->effectiveTransportCoefficient = (3.0 * absorptionCoefficient * reducedExtinctionCoefficient ).Sqrt();	// sigma_tr = sqrt(3 * sigma_a * sigma_t')
	// float F = Ft(relativeRefractiveIndex, wo) Ft(relativeRefractiveIndex, wi)
}


FluorescentMaterial::~FluorescentMaterial() {
}


// Fluorescent Method Definitions
BSDF *FluorescentMaterial::GetBSDF(const DifferentialGeometry &dgGeom,
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
	string fluoroFile = "white.txt";
	float relRefractiveIndex = 1.3f;
	float defaultSigmaSPrime[3] = {2.29f, 2.39f, 1.97f};
	float defaultSigmaA[3] = {0.003f, 0.0034f, 0.046f};
	Spectrum reducedScatteringCoefficient = Spectrum(defaultSigmaSPrime);
	Spectrum absorptionCoefficient = Spectrum(defaultSigmaA);
	reducedScatteringCoefficient.setValueAtWavelength(5.0f, 370);
	reducedScatteringCoefficient.setValueAtWavelength(4.0f, 380);
	reducedScatteringCoefficient.setValueAtWavelength(4.0f, 360);
	absorptionCoefficient.setValueAtWavelength(0.003f, 370);
	absorptionCoefficient.setValueAtWavelength(0.004f, 380);
	absorptionCoefficient.setValueAtWavelength(0.004f, 360);
	if(paramSet != NULL)
	{
		printf("\tUSING PARAM SET\n");
		fluoroFile = paramSet->FindOneString("reradiation", fluoroFile);
		reducedScatteringCoefficient = paramSet->FindOneSpectrum("reducedscatteringcoefficient", reducedScatteringCoefficient);
		absorptionCoefficient = paramSet->FindOneSpectrum("absorptioncoefficient", absorptionCoefficient);
		relRefractiveIndex= paramSet->FindOneFloat("relativerefractiveindex", relRefractiveIndex);
	}
	
	
	return new FluorescentMaterial(Kd, Ks, roughness, bumpMap, fluoroFile, reducedScatteringCoefficient, absorptionCoefficient, relRefractiveIndex);
}
