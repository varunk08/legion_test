//-------------------------------------------------------------------------------
///
/// \file       materials.h
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    1.0
/// \date       September 2, 2013
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------

#ifndef _MATERIALS_H_INCLUDED_
#define _MATERIALS_H_INCLUDED_

#include "scene.h"

//-------------------------------------------------------------------------------

class MtlBlinn : public Material
{

 public:
 MtlBlinn() : diffuse(0.5f,0.5f,0.5f), specular(0.7f,0.7f,0.7f), glossiness(20.0f),
    reflection(0,0,0), refraction(0,0,0), absorption(0,0,0), ior(1) {}

  void SetDiffuse(cyColor dif) { diffuse = dif; }
  void SetSpecular(cyColor spec) { specular = spec; }
  void SetGlossiness(float gloss) { glossiness = gloss; }
  void SetReflection(cyColor reflect) { reflection = reflect; }
  void SetRefraction(cyColor refract) { refraction = refract; }
  void SetAbsorption(cyColor absorp ) { absorption = absorp; }
  void SetRefractionIndex(float _ior) { ior = _ior; }

  cyColor GetDiffuse() {return  diffuse;}
  cyColor GetSpecular() {return  specular; }
  float GetGlossiness() {return  glossiness; }
  cyColor GetReflection() {return  reflection; }
  cyColor GetRefraction() {return  refraction; }
  cyColor GetAbsorption() {return  absorption; }
  float GetRefractionIndex() {return  ior; }

  //Main shade function to be filled in later
  cyColor Shade(const Ray &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount) const
    {
      cyColor shade;
      cyColor rShade = cyColor(0,0,0);
      cyColor tShade = cyColor(0,0,0);
      const Material *mat;
      mat = hInfo.node->GetMaterial();
      const MtlBlinn* mb =static_cast<const MtlBlinn*>(mat);

      Point3 P;
      P.Set(hInfo.p.x,hInfo.p.y,hInfo.p.z);
      Ray iRay = ray;
      cyColor ambInt = mb->diffuse;
      cyColor allOther = cyColor(0,0,0);
      cyColor diffuse = mb->diffuse;;
      cyColor ambComponent = cyColor(0,0,0);
    
      for ( unsigned int i=0; i<lights.size(); i++ ) {
        if(lights[i]->IsAmbient()){
	  cyColor intensity = lights[i]->Illuminate(hInfo.p);
	  ambComponent += (ambInt * intensity);
	  continue;
        }
        else{
	  Point3 L = -lights[i]->Direction(P);
	  L.Normalize();

	  Point3 V = ray.p - P;
	  V.Normalize();

	  Point3 LplusV = L + V;
	  Point3 H = (L+V)/LplusV.Length();
	  H.Normalize();
            
	  float alpha = mb->glossiness;
	  Point3 N = hInfo.N;
	  float S = H.Dot(N);
	  S = pow((float)S,alpha);
	  float costheta = L.Dot(N)/(L.Length() * N.Length());
	  cyColor intensity = lights[i]->Illuminate(P);

	  allOther += intensity * (costheta>0?costheta:0) * (diffuse + S * (mb->specular)) ;
        }
        /* finally add inta*cola + intall*costheta*(cold + s* colS)*/
        shade = ambComponent  + allOther;
      }

      return shade;

    }    
  ~MtlBlinn(){};

    
 private:
  cyColor diffuse, specular, reflection, refraction;
  float glossiness;
  cyColor absorption;
  float ior;	// index of refraction
};

//-------------------------------------------------------------------------------

#endif
