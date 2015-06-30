//
//  Plane.h
//  RayTracePrj5
//
//  Created by Varun kumar Karuppannan on 28/09/13.
//  Copyright (c) 2013 Varun kumar Karuppannan. All rights reserved.
//

#include <iostream>
#include <string.h>
#include <strings.h>
#include <cmath>
#include "scene.h"
#include "cyTriMesh.h"
#include "cyBVH.h"
using namespace std;

class TriObj : public Object, private cyTriMesh
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
    void SetTransform(const Matrix3 &nodeToWorld, const Matrix3 &itm,  const Point3 &pos){}   ;

    bool IntersectRay( const Ray &ray, HitInfo& hInfo, int hitSide=HIT_FRONT )const
    {
        bool hit = false;
        if(TraceBVHNode(ray, hInfo, hitSide, bvh.GetRootNodeID())){
            hit = true;
       /*
        //get bound box
        Box aabb = GetBoundBox();
        
        //check interseciton with bound box
        bool boxHit = aabb.IntersectRay(ray, BIGFLOAT);
        //then check for interseciton with each triangle - already implemented
        if(boxHit){
            Point3 P(ray.p.x,ray.p.y,ray.p.z);
            Point3 D(ray.dir.x,ray.dir.y,ray.dir.z);
            Point3 N;
            
            for(unsigned int i = 0; i < NF(); i++){
                if(IntersectTriangle(ray, hInfo, hitSide, i)){
                    hit = true;
                    
                }
            }
        }*/
        }
        return hit;
        
    }
    
    Box GetBoundBox() const { return Box(GetBoundMin(),GetBoundMax()); }
    
	bool Load(const char *filename)
	{
		if ( ! LoadFromFileObj( filename ) ) return false;
		if ( ! HasNormals() ) ComputeNormals();
		ComputeBoundingBox();
        bvh.SetMesh(this,4);
		return true;

    }
private:
    cyBVHTriMesh bvh;
    bool TraceBVHNode( const Ray &ray, HitInfo &hInfo, int hitSide, unsigned int nodeID ) const
    {
        
        Box aabb = Box(bvh.GetNodeBounds(nodeID));
        bool triangleHit = false;
        if(!aabb.IntersectRay(ray, BIGFLOAT)){
            return false;
        }
        else{
           
            
            
            if(bvh.IsLeafNode(nodeID)){

                unsigned int* elements = (unsigned int*)bvh.GetNodeElements(nodeID);
                for(unsigned int i = 0; i < bvh.GetNodeElementCount(nodeID); i++){
                    if(IntersectTriangle(ray, hInfo, hitSide, elements[i])){
                        triangleHit = true;
                    }
                }
                return triangleHit;
            }
            
            else{
                if(TraceBVHNode(ray, hInfo, hitSide, bvh.GetFirstChildNode(nodeID))){
                    triangleHit = true;
                }
                if(TraceBVHNode(ray, hInfo, hitSide, bvh.GetSecondChildNode(nodeID))){
                    triangleHit = true;
                }
            }
            
        }
        return triangleHit;

    }
    
    bool IntersectTriangle( const Ray &ray, HitInfo &hInfo, int hitSide, unsigned int faceID ) const
    {
        //get the vertices of each face with the faceID
        Point3 A, B, C, P, D, N, H;
        float t, a, a1, a2, alpha, beta;
        float bias = 1e-6f;
        Point2 A2d, B2d, C2d, H2d;
        A = V(F(faceID).v[0]);
        B = V(F(faceID).v[1]);
        C = V(F(faceID).v[2]);
	//if(Box_IntersectTriangle(ray, hInfo, A, B, C)){
	// return true;
	//}

        //Do Ray triangle intersection with the three vertices thus obtained
        P = ray.p;
        D = ray.dir;
        Point3 nTemp = (B - A).Cross(C - A);
        //N = nTemp/ nTemp.Length();
        nTemp.Normalize();
        N = nTemp;

        if(N.Dot(P-A) < bias) return false; /* we are in the back face */
        if(N.Dot(D) == 0) return false; /* ray is parallel */
        t = N.Dot(P-A)/-(N.Dot(D));
        
        //if the t-value is less than what is there in hit info then assign and signal hit
        if(t < hInfo.z && t >= bias && t< BIGFLOAT){
            H = P + t * D; /* point of intersection on the plane */
            int projAxis = 0; //0-xaxis, 1-yaxis, 2-zaxis
            //project H onto an axis aligned plane - find the maximum component
            if(abs(H.x) > abs(H.y)){
                if(abs(H.x) > abs(H.z)){
                    projAxis = 0;
                }
                else{
                    projAxis = 2;
                }
            }
            else if(abs(H.y) > abs(H.z)){
                projAxis = 1;
            }
            else{
                projAxis = 2;
            }
            //cout<<projAxis<<endl;
            switch (projAxis) {
                case 0: /* project onto the x-axis*/
                    A2d = Point2(A.y, A.z);
                    B2d = Point2(B.y, B.z);
                    C2d = Point2(C.y, C.z);
                    H2d = Point2(H.y, H.z);
                    break;
                case 1:/* project onto the y-axis*/
                    A2d = Point2(A.x, A.z);
                    B2d = Point2(B.x, B.z);
                    C2d = Point2(C.x, C.z);
                    H2d = Point2(H.x, H.z);
                    break;
                case 2:/* project onto the z-axis*/
                    A2d = Point2(A.x, A.y);
                    B2d = Point2(B.x, B.y);
                    C2d = Point2(C.x, C.y);
                    H2d = Point2(H.x, H.y);
                default:
                    break;
            }
            
            a = (B2d - A2d).Cross(C2d - A2d)/ 2.0;
            a1 = (B2d - A2d).Cross(H2d - A2d)/ 2.0;
            a2 = (H2d - A2d).Cross(C2d - A2d)/ 2.0;
            alpha = a1/a;
            beta = a2/a;
            if(alpha < -bias || beta < -bias || alpha + beta > 1+bias){
                return false; /* point H is outside triangle */
            }
            float b0 = 1 - beta - alpha;
            float b1 = beta;
            float b2 = alpha;
            H = b0 * A + b1 * B + b2 * C;
            Point3 Normal = GetNormal(faceID, Point3(b0,b1,b2));
            N = Normal;

            hInfo.front = true;
            hInfo.p = H;
            hInfo.N = N;
            hInfo.z = t;
            return true;
            
        }
        
        return false;
    }

  bool Box_IntersectTriangle(const Ray &ray, HitInfo &hInfo, Point3 p1, Point3 p2, Point3 p3) const
  {
    Point3 A, B, C, N, P, D, H;
    float bias = 1e-6f;
    float t, a, a1, a2, alpha, beta;
    Point2 A2d, B2d, C2d, H2d;

    A = p1; B = p2; C = p3;
    P = ray.p;
    D = ray.dir;
    Point3 nTemp = (B-A).Cross(C-A);
    nTemp.Normalize();
    N = nTemp;
    
    if(N.Dot(P-A) < bias) return false;
    if(N.Dot(D) == 0) return false;
    t = N.Dot(P-A)/ -(N.Dot(D));

  if( t < hInfo.z && t >= bias && t < BIGFLOAT ) 
    {
      H = P + t * D;
      int projAxis = 0;
            if(abs(H.x) > abs(H.y)){
                if(abs(H.x) > abs(H.z)){
                    projAxis = 0;
                }
                else{
                    projAxis = 2;
                }
            }
            else if(abs(H.y) > abs(H.z)){
                projAxis = 1;
            }
            else{
                projAxis = 2;
            }
            //cout<<projAxis<<endl;
            switch (projAxis) {
                case 0: /* project onto the x-axis*/
                    A2d = Point2(A.y, A.z);
                    B2d = Point2(B.y, B.z);
                    C2d = Point2(C.y, C.z);
                    H2d = Point2(H.y, H.z);
                    break;
                case 1:/* project onto the y-axis*/
                    A2d = Point2(A.x, A.z);
                    B2d = Point2(B.x, B.z);
                    C2d = Point2(C.x, C.z);
                    H2d = Point2(H.x, H.z);
                    break;
                case 2:/* project onto the z-axis*/
                    A2d = Point2(A.x, A.y);
                    B2d = Point2(B.x, B.y);
                    C2d = Point2(C.x, C.y);
                    H2d = Point2(H.x, H.y);
                default:
                    break;
            }
            
            a = (B2d - A2d).Cross(C2d - A2d)/ 2.0;
            a1 = (B2d - A2d).Cross(H2d - A2d)/ 2.0;
            a2 = (H2d - A2d).Cross(C2d - A2d)/ 2.0;
            alpha = a1/a;
            beta = a2/a;
            if(alpha < -bias || beta < -bias || alpha + beta > 1+bias){
                return false; /* point H is outside triangle */
            }
            float b0 = 1 - beta - alpha;
            float b1 = beta;
            float b2 = alpha;
            H = b0 * A + b1 * B + b2 * C;

	    hInfo.front = true;
	    hInfo.p = H;
	    hInfo.N = N;
	    hInfo.z = t;

	    return true;
	    
    }
  return false;
    
  }

};
