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

  //Main shade function to be filled in later
  cyColor Shade(const Ray &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount) const
    {
      return cyColor(1,1,1);
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
