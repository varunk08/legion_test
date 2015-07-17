/* Volume data calculations using Legion */
#ifndef VOLUME_GEOM_H
#define VOLUME_GEOM_H

//-----------------------------------------------------------------------------------------------------

#include <vector>
#include <iostream>
#include <string.h>
#include <cmath>
#include "scene.h"
#include "common_types.h"
#include "cyColor.h"
#include "utils.h"
#include "cyMatrix3.h"
#include "lights.h"
#include "legion.h"

//-----------------------------------------------------------------------------------------------------
using namespace std;
using namespace LegionRuntime::HighLevel;
using namespace LegionRuntime::Accessor;

//-----------------------------------------------------------------------------------------------------

class VolumeGeometry : public Object
{


 private:
  Point3 pmin;
  Point3 pmax;
  int tf_size;
  uchar minData;
  uchar maxData;
  LightList lights;
  Matrix3 tm, itm;
  Point3 worldPos;
  std::vector<PhysicalRegion> vol_regions;
  Context ctx;
  HighLevelRuntime *runtime;
  int data_xdim;
  int data_ydim;
  int data_zdim;
  GPoint3fRA acc_corner_pos;
  GucharRA acc_vol_data;
  GPoint3fRA acc_gradients;
  GRGBColorRA acc_color_tf;
  GfloatRA acc_alpha_tf;
  LogicalRegion corner_pos_lr;
  LogicalRegion tf_lr;


 public:
  //constructor
  VolumeGeometry(int x_width, int y_width, int z_width,
		 const std::vector<PhysicalRegion> &regions, Context ctx, HighLevelRuntime *runtime);


  bool IntersectRay(const Ray &ray, HitInfo &hInfo, int hitSide=HIT_FRONT ) const;
  Box GetBoundBox() const;
  void init_logical_regions(LogicalRegion vol_data_lr, LogicalRegion tf_lr,
			    GucharRA &vol_acc, GPoint3fRA &acc, GPoint3fRA &acc_grads, GRGBColorRA &acc_col,
			    GfloatRA &acc_alpha);
  void init_tf_bounds(int min, int max);
  void set_lights(LightList &lList);
  void SetTransform(const Matrix3 &nodeToWorld, const Matrix3 &n2w_itm, const Point3 &pos );

 private:

  bool Box_IntersectTriangle(const Ray &ray, HitInfo &hInfo, Point3 p1, Point3 p2, Point3 p3, 
			     bool raymarch=false) const;
  Point3 EstimateGradient(int x, int y, int z) const;
  float TrilinearInterpolate(Point3 samplePt, int x, int y, int z) const;
  inline int getIndex(int x, int y, int z) const;
  float GetAlpha(float density) const;
  cyColor GetColor(float density) const;
  cyColor DoPhongShading(const Ray &ray,const Point3 &pt, const Point3 &norm, cyColor preCol) const;
  cyColor RayMarch(const Ray& ray, HitInfo &hInfo, HitInfo& start, HitInfo& end, bool &noData) const;
  Point3 TransposeMult( const Matrix3 &m, const Point3 &dir ) const;
  Point3 TransformFrom( const Point3 &p ) const;
  Point3 VectorTransformFrom( const Point3 &dir ) const;


};
//-----------------------------------------------------------------------------------------------------
#endif
