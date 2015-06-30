#ifndef XMLLOAD_H
#define XMLLOAD_H

#include "scene.h"
#include "objects.h"
#include "sphere.h"
#include "TriObj.h"
#include "plane.h"
#include "materials.h"
#include "lights.h"
#include "tinyxml/tinyxml.h"
#include "boxobject.h"

//int LoadScene(const char *filename);
int LoadScene(const char *filename,
	      Node         &rootNode,
	      Camera       &camera,
	      RenderImage  &renderImage,
	      MaterialList &materials,
	      LightList    &lights,
	      ObjFileList  &objList,
	      Sphere       &theSphere,
	      BoxObject    &theBoxObject,
	      Plane        &thePlane);

void LoadScene(TiXmlElement *element,
	       Node &rootNode,
	       ObjFileList &objList,
	       MaterialList &materials,
	       LightList &lights,
	       Sphere       &theSphere,
	       BoxObject    &theBoxObject,
	       Plane        &thePlane);
void LoadNode(Node *node,
	      TiXmlElement *element,
	      ObjFileList &objList,
	      Sphere       &theSphere,
	      BoxObject    &theBoxObject,
	      Plane        &thePlane,
	      int level=0);
void LoadTransform( Transformation *trans, TiXmlElement *element, int level );
void LoadMaterial(TiXmlElement *element, MaterialList &materials);
void LoadLight(TiXmlElement *element, LightList &lights);
void ReadVector(TiXmlElement *element, Point3 &v);
void ReadColor (TiXmlElement *element, cyColor  &c);
void ReadFloat (TiXmlElement *element, float  &f, const char *name="value");


#endif
