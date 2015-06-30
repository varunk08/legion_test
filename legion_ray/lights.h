//-------------------------------------------------------------------------------
///
/// \file       lights.h
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    2.0
/// \date       September 9, 2013
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------

#ifndef _LIGHTS_H_INCLUDED_
#define _LIGHTS_H_INCLUDED_

#include "scene.h"

//-------------------------------------------------------------------------------

class GenLight : public Light
{
protected:
  /*	void SetViewportParam(int lightID, cyColorA ambient, cyColorA intensity, Point4 pos ) const{
        
        glEnable ( GL_LIGHT0 + lightID );
        glLightfv( GL_LIGHT0 + lightID, GL_AMBIENT,  &ambient.r );
        glLightfv( GL_LIGHT0 + lightID, GL_DIFFUSE,  &intensity.r );
        glLightfv( GL_LIGHT0 + lightID, GL_SPECULAR, &intensity.r );
        glLightfv( GL_LIGHT0 + lightID, GL_POSITION, &pos.x );
        
	};*/
  static float Shadow(Ray ray, float t_max=BIGFLOAT);
};

//-------------------------------------------------------------------------------

class AmbientLight : public GenLight
{
public:
	AmbientLight() : intensity(0,0,0) {}
	virtual cyColor Illuminate(Point3 p) const { return intensity; }
	virtual Point3 Direction(Point3 p) const { return Point3(0,0,0); }
	virtual bool IsAmbient() const { return true; }
	//	virtual void SetViewportLight(int lightID) const { SetViewportParam(lightID,intensity,0.0f,Point4(0,0,0,1)); }
    
	void SetIntensity(cyColor intens) { intensity=intens; }
private:
	cyColor intensity;
};

//-------------------------------------------------------------------------------

class DirectLight : public GenLight
{
public:
	DirectLight() : intensity(0,0,0), direction(0,0,1) {}
	virtual cyColor Illuminate(Point3 p) const { return Shadow(Ray(p,-direction)) * intensity; }
	virtual Point3 Direction(Point3 p) const { return direction; }
	//	virtual void SetViewportLight(int lightID) const { SetViewportParam(lightID,0.0f,intensity,Point4(-direction,0.0f)); }
    
	void SetIntensity(cyColor intens) { intensity=intens; }
	void SetDirection(Point3 dir) { direction=dir.GetNormalized(); }
private:
	cyColor intensity;
	Point3 direction;
};

//-------------------------------------------------------------------------------

class PointLight : public GenLight
{
public:
	PointLight() : intensity(0,0,0), position(0,0,0) {}
	virtual cyColor Illuminate(Point3 p) const { return Shadow(Ray(p,position-p),1) * intensity; }
	virtual Point3 Direction(Point3 p) const { return (p-position).GetNormalized(); }
	//	virtual void SetViewportLight(int lightID) const { SetViewportParam(lightID,0.0f,intensity,Point4(position,1.0f)); }
	void SetIntensity(cyColor intens) { intensity=intens; }
	void SetPosition(Point3 pos) { position=pos; }
    
private:
	cyColor intensity;
	Point3 position;
};

//-------------------------------------------------------------------------------

#endif
