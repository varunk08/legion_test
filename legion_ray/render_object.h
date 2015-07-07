#ifndef _RENDER_OBJECT
#define _RENDER_OBJECT
#include <vector>
#include "scene.h"
#include "materials.h"
#include "sphere.h"
class RenderObject: public Transformation
{
 public:
  Sphere theObject;
  MtlBlinn theMaterial;
  std::vector<Light> theLightList;

  //creates a copy of the object
  void SetObject (const Sphere &object) { this->theObject = object; }
  void SetMaterial (const Material* mtl) {
    //this->theMaterial = new MtlBlinn();
    this->theMaterial.SetDiffuse(((MtlBlinn*)mtl)->GetDiffuse()); 
    this->theMaterial.SetSpecular(((MtlBlinn*)mtl)->GetSpecular());
    this->theMaterial.SetGlossiness(((MtlBlinn*)mtl)->GetGlossiness());
    this->theMaterial.SetReflection(((MtlBlinn*)mtl)->GetReflection());
    this->theMaterial.SetRefraction(((MtlBlinn*)mtl)->GetRefraction());
    this->theMaterial.SetAbsorption(((MtlBlinn*)mtl)->GetAbsorption());
    this->theMaterial.SetRefractionIndex(((MtlBlinn*)mtl)->GetRefractionIndex());
  }
  void SetLightList(LightList lightList) {
    this->theLightList.clear();
    for(int i=0; i< lightList.size(); i++){
      if(lightList.at(i)){
	this->theLightList.at(i) = *(lightList.at(i));
      }
    }

  }
  void SetTransformation(Matrix3 m) { this->Transform(m); };
  // Transformations
  Ray ToNodeCoords( const Ray &ray ) const
  {
    Ray r;
    r.p   = TransformTo(ray.p);
    r.dir = TransformTo(ray.p + ray.dir) - r.p;
    return r;
  }
  void FromNodeCoords( HitInfo &hInfo ) const
  {
    hInfo.p = TransformFrom(hInfo.p);
    hInfo.N = VectorTransformFrom(hInfo.N).GetNormalized();
  }


};


#endif
