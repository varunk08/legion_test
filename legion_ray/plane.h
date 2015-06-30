//
//  Plane.h
//  RayTracePrj5
//
//  Created by Varun kumar Karuppannan on 28/09/13.
//  Copyright (c) 2013 Varun kumar Karuppannan. All rights reserved.
//
#include "scene.h"
#include <iostream>
#include <string.h>
#include <strings.h>
#include <cmath>
//#include <GLUT/GLUT.h>

using namespace std;

class Plane : public Object
{
    
    
private:
    char* name;
public:
    
    void setName(const char* newName){
        cout<<newName<<endl;
        if ( newName ) {
            int n = strlen(newName);
            name = new char[n+1];
            for ( int i=0; i<n; i++ ) name[i] = newName[i];
            name[n] = '\0';
        } else { name = NULL; }
    };
    Box GetBoundBox() const { return Box(-1,-1,0,1,1,0); }
    void SetTransform(const Matrix3 &nodeToWorld,const Matrix3 &itm, const Point3 &pos) {}
    bool IntersectRay( const Ray &ray, HitInfo& hInfo, int hitSide=HIT_FRONT )const
    {
        
        float bias = 1e-6f;
        Point3 P(ray.p.x,ray.p.y,ray.p.z);
        Point3 D(ray.dir.x,ray.dir.y,ray.dir.z);
        Point3 N(0,0,1);
        
        float t = -(P.z/D.z);
        if(t >= 0 && t < BIGFLOAT && t<hInfo.z){
            Point3 pHit(P.x + t * D.x, P.y + t * D.y, P.z + t * D.z);
            if(pHit.x >= -1 && pHit.x <= 1 && pHit.y >=-1 && pHit.y <=1){
                hInfo.z = t;
//                cout<<"pHit: "<<pHit.x<<" pHit: "<<pHit.y<<" pHit: "<<pHit.z<<endl;
                hInfo.p.Set(pHit.x, pHit.y, pHit.z);
                hInfo.N = N;
                if(N.Dot(D) < 0){
                    hInfo.front = false;
                }
                else hInfo.front = true;
                return true;
            }
        }
        
//        
//        Point3 Q = Point3(1,1,0);
//        Point3 P = ray.p;
//        Point3 D = ray.dir;
//        Point3 N = Point3(0,0,1);
//        
//        float t = N.Dot(P - Q) / -(N.Dot(D));
//        
//        if(t >= bias &&t < BIGFLOAT){
//                    hInfo.z = t;
//                    hInfo.p = Point3(P.x + t * D.x, P.y + t * D.y, 0);
//                    hInfo.N = Point3(hInfo.p.x,hInfo.p.y,1);
//                    hInfo.front = true;
//                    return true;
//        }
        return false;
        
    }
    
    /*    void ViewportDisplay() const
    {
        const int resolution = 32;
        glPushMatrix();
        glScalef(2.0f/resolution,2.0f/resolution,2.0f/resolution);
        glNormal3f(0,0,1);
        glBegin(GL_QUADS);
        for ( int y=0; y<resolution; y++ ) {
            float yy = y - resolution/2;
            for ( int x=0; x<resolution; x++ ) {
                float xx = x - resolution/2;
                glVertex3f( yy,   xx,   0 );
                glVertex3f( yy+1, xx,   0 );
                glVertex3f( yy+1, xx+1, 0 );
                glVertex3f( yy,   xx+1, 0 );
            }
        }
        glEnd();
        glPopMatrix();
    }
    
    */
};

