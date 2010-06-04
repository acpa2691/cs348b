//Modified subsurface integrator from [song,yu]
//Outfitted with sophisticated spectral scattering effects
// subsurface.cpp*
#include "pbrt.h"
#include "color.h"
#include "transport.h"
#include "scene.h"
#include "mc.h"
#include "primitive.h"
#include "octree.h"
#include "material.h"
#include "../materials/fluorescent.cpp"
#include "../shapes/trianglemesh.cpp"

// DirectLighting Declarations
enum LightStrategy { SAMPLE_ALL_UNIFORM, SAMPLE_ONE_UNIFORM,
	SAMPLE_ONE_WEIGHTED  // NOBOOK
};

class SubsurfaceIntegrator : public SurfaceIntegrator {
public:

	SubsurfaceIntegrator()
	{
	}
	
	SubsurfaceIntegrator(LightStrategy ls, int md)
	{
		printf("Making subsurface integrator w00t\n");
		maxDepth = md;
		rayDepth = 0;
		strategy = ls;
		avgY = avgYsample = cdf = NULL;
		overallAvgY = 0.;
	}
	
	~SubsurfaceIntegrator();
	Spectrum Li(const Scene *scene, const RayDifferential &ray, const Sample *sample, float *alpha) const;
	void RequestSamples(Sample *sample, const Scene *scene);
	void Preprocess(const Scene * scene);
private:
	// DirectLighting Private Data
	LightStrategy strategy;
	mutable int rayDepth; // NOBOOK
	int maxDepth; // NOBOOK
	// Declare sample parameters for light source sampling
	int *lightSampleOffset, lightNumOffset;
	int *bsdfSampleOffset, *bsdfComponentOffset;
	mutable float *avgY, *avgYsample, *cdf;
	mutable float overallAvgY;
	
	
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

SubsurfaceIntegrator::~SubsurfaceIntegrator()
{
	delete[] avgY;
	delete[] avgYsample;
	delete[] cdf;
}

Spectrum SubsurfaceIntegrator::Li(const Scene *scene, const RayDifferential &ray, const Sample *sample, float *alpha) const
{
	/*
	Pseudocode
		
	if ray intersects scene
		L = different components
		for now, assume just subsurface from direct lighting component
		
		Compute subsurface:
			
			get BSSRDF
			integrate bssrdf over area
			for now, assume bssrdf constant over area (uniform lighting)
			
			get material properties
			have scene, intersection, primitive
			
			if (! isect->p->hasBSSRDF / isFluorescent / hasSubsurfaceScattering):
				use regular brdf stuff
				
			isect->getBSSRDF
			
			uniform sample all lights
	*/

	Intersection isect;
	Spectrum L(0.);
	if (scene->Intersect(ray, &isect)) {
		if (alpha) *alpha = 1.;
		// Evaluate BSDF at hit point
		BSDF *bsdf = isect.GetBSDF(ray);
		Vector wo = -ray.d;
		const Point &p = bsdf->dgShading.p;
		const Normal &n = bsdf->dgShading.nn;
		// Compute emitted light if ray hit an area light source
		L += isect.Le(wo);
		// Compute direct lighting for _DirectLighting_ integrator
		if (scene->lights.size() > 0) {
			// Apply subsurface scattering
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
						if (dr < 0.8*lu) continue;
						
						Point Pv = xi + lu * (1 + 4*A/3)* n;
						float dv = (xo - Pv).Length();
						Spectrum C1 = lu * (sigmatr + 1 / dr);
						Spectrum C2 = zv * (sigmatr + 1 / dv);
						Spectrum diffuseReflectance = (0.25 * INV_PI) * alphaprime * (C1 * Exp( -dr*sigmatr ) / (dr*dr) + C2 * Exp( -dv*sigmatr ) / (dv*dv));
						//printf("i have diffuse reflectance ");
						//diffuseReflectance.printSelf();
						//Spectrum diffuseReflectance = (0.25 * INV_PI) * alphaprime * (C1 * Exp( -dr*sigmatr ) + C2 * Exp( -dv*sigmatr ));
						//Spectrum diffuseReflectance = (0.25 * INV_PI) * alphaprime * (C1 * Exp( -dr*sigmatr ) / dr + C2 * Exp( -dv*sigmatr ) / dv);
						
						Spectrum out = fluoro->output(E);
						//printf("adding irradiance");
						//out.printSelf();
						L += diffuseReflectance * Area * out; // * Ft
					}
				}
				L *= AbsDot(wo, n);
			}
			
			if (bsdf->NumComponents(BxDFType(BSDF_REFLECTION)) == 1) {
				L = L*0.6 + 0.4*Spectrum(0.9f);
			}
			
			// Apply direct lighting strategy
			if(!isect.primitive->IsFluorescent()) {
				switch (strategy) {
						case SAMPLE_ALL_UNIFORM:
							L += UniformSampleAllLights(scene, p, n, wo, bsdf,
								sample, lightSampleOffset, bsdfSampleOffset,
								bsdfComponentOffset);
							break;
						case SAMPLE_ONE_UNIFORM:
							L += UniformSampleOneLight(scene, p, n, wo, bsdf,
								sample, lightSampleOffset[0], lightNumOffset,
								bsdfSampleOffset[0], bsdfComponentOffset[0]);
							break;
						case SAMPLE_ONE_WEIGHTED:
							L += WeightedSampleOneLight(scene, p, n, wo, bsdf,
								sample, lightSampleOffset[0], lightNumOffset,
								bsdfSampleOffset[0], bsdfComponentOffset[0], avgY,
								avgYsample, cdf, overallAvgY);
							break;
					}
			}
			//*/
			
		}
		
		return L;
		if (rayDepth++ < maxDepth) {
			//printf("at depth: %d\n", rayDepth);
			Vector wi;
			// Trace rays for specular reflection and refraction
			Spectrum f = bsdf->Sample_f(wo, &wi,
				BxDFType(BSDF_REFLECTION | BSDF_SPECULAR));
			if (!f.Black()) {
				// Compute ray differential _rd_ for specular reflection
				RayDifferential rd(p, wi);
				rd.hasDifferentials = true;
				rd.rx.o = p + isect.dg.dpdx;
				rd.ry.o = p + isect.dg.dpdy;
				// Compute differential reflected directions
				Normal dndx = bsdf->dgShading.dndu * bsdf->dgShading.dudx +
					bsdf->dgShading.dndv * bsdf->dgShading.dvdx;
				Normal dndy = bsdf->dgShading.dndu * bsdf->dgShading.dudy +
					bsdf->dgShading.dndv * bsdf->dgShading.dvdy;
				Vector dwodx = -ray.rx.d - wo, dwody = -ray.ry.d - wo;
				float dDNdx = Dot(dwodx, n) + Dot(wo, dndx);
				float dDNdy = Dot(dwody, n) + Dot(wo, dndy);
				rd.rx.d = wi -
				          dwodx + 2 * Vector(Dot(wo, n) * dndx +
						  dDNdx * n);
				rd.ry.d = wi -
				          dwody + 2 * Vector(Dot(wo, n) * dndy +
						  dDNdy * n);
				L += scene->Li(rd, sample) * f * AbsDot(wi, n);
			}
			f = bsdf->Sample_f(wo, &wi,
				BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR));
			//printf("got sample\n");
			if (!f.Black()) {
				// Compute ray differential _rd_ for specular transmission
				RayDifferential rd(p, wi);
				rd.hasDifferentials = true;
				rd.rx.o = p + isect.dg.dpdx;
				rd.ry.o = p + isect.dg.dpdy;
				
				float eta = bsdf->eta;
				Vector w = -wo;
				if (Dot(wo, n) < 0) eta = 1.f / eta;
				
				Normal dndx = bsdf->dgShading.dndu * bsdf->dgShading.dudx + bsdf->dgShading.dndv * bsdf->dgShading.dvdx;
				Normal dndy = bsdf->dgShading.dndu * bsdf->dgShading.dudy + bsdf->dgShading.dndv * bsdf->dgShading.dvdy;
				
				Vector dwodx = -ray.rx.d - wo, dwody = -ray.ry.d - wo;
				float dDNdx = Dot(dwodx, n) + Dot(wo, dndx);
				float dDNdy = Dot(dwody, n) + Dot(wo, dndy);
				
				float mu = eta * Dot(w, n) - Dot(wi, n);
				float dmudx = (eta - (eta*eta*Dot(w,n))/Dot(wi, n)) * dDNdx;
				float dmudy = (eta - (eta*eta*Dot(w,n))/Dot(wi, n)) * dDNdy;
				
				rd.rx.d = wi + eta * dwodx - Vector(mu * dndx + dmudx * n);
				rd.ry.d = wi + eta * dwody - Vector(mu * dndy + dmudy * n);
				L += scene->Li(rd, sample) * f * AbsDot(wi, n);
			}
		}
		--rayDepth;
	}
	else {
		// Handle ray with no intersection
		if (alpha) *alpha = 0.;
		for (u_int i = 0; i < scene->lights.size(); ++i)
			L += scene->lights[i]->Le(ray);
		if (alpha && !L.Black()) *alpha = 1.;
		return L;
	}
	return L;
}

void SubsurfaceIntegrator::RequestSamples(Sample *sample, const Scene *scene)
{
	if (strategy == SAMPLE_ALL_UNIFORM) {
		// Allocate and request samples for sampling all lights
		u_int nLights = scene->lights.size();
		lightSampleOffset = new int[nLights];
		bsdfSampleOffset = new int[nLights];
		bsdfComponentOffset = new int[nLights];
		for (u_int i = 0; i < nLights; ++i) {
			const Light *light = scene->lights[i];
			int lightSamples =
				scene->sampler->RoundSize(light->nSamples);
			lightSampleOffset[i] = sample->Add2D(lightSamples);
			bsdfSampleOffset[i] = sample->Add2D(lightSamples);
			bsdfComponentOffset[i] = sample->Add1D(lightSamples);
		}
		lightNumOffset = -1;
	}
	else {
		// Allocate and request samples for sampling one light
		lightSampleOffset = new int[1];
		lightSampleOffset[0] = sample->Add2D(1);
		lightNumOffset = sample->Add1D(1);
		bsdfSampleOffset = new int[1];
		bsdfSampleOffset[0] = sample->Add2D(1);
		bsdfComponentOffset = new int[1];
		bsdfComponentOffset[0] = sample->Add1D(1);
	}
}

void SubsurfaceIntegrator::Preprocess(const Scene * scene)
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
			cacheEntry.E = EstimateIrradiance(scene, cacheEntry.p, cacheEntry.n, lightSampleOffset);
		}
	}
}

extern "C" DLLEXPORT SurfaceIntegrator *CreateSurfaceIntegrator(const ParamSet &params) {
	int maxDepth = params.FindOneInt("maxdepth", 5);
	LightStrategy strategy;
	string st = params.FindOneString("strategy", "all");
	if (st == "one") strategy = SAMPLE_ONE_UNIFORM;
	else if (st == "all") strategy = SAMPLE_ALL_UNIFORM;
	else if (st == "weighted") strategy = SAMPLE_ONE_WEIGHTED;
	else {
		Warning("Strategy \"%s\" for direct lighting unknown. "
			"Using \"all\".", st.c_str());
		strategy = SAMPLE_ALL_UNIFORM;
	}
	return new SubsurfaceIntegrator(strategy, maxDepth);
}
