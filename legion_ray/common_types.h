#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H
//-----------------------------------------------------------------------------------------------------
#include "legion.h"
#include "scene.h"
//-----------------------------------------------------------------------------------------------------
enum TaskID{
  TOP_LEVEL_ID,
  PER_PIXEL_TASK_ID,
  CHECK_TASK_ID,
  VOLUME_RENDER_TASK_ID,
};
enum FieldSpaceID{
  FID_OUT,
  FID_VOL_DATA,
  FID_CORNER_POS,
  FID_GRADIENTS,
  FID_COLOR_TF,
  FID_ALPHA_TF,
};

//Info needed by each task to perform rendering calculations
struct SceneData 
{
  Ray ray; //ray info
  int coords[2];
  int index;
};

struct RGBColor
{
  float r; float g; float b;
};

struct Point3f
{
  float x; float y; float z;
};

typedef LegionRuntime::Accessor::RegionAccessor < LegionRuntime::Accessor::AccessorType::Generic, Point3f> GPoint3fRA;
typedef LegionRuntime::Accessor::RegionAccessor < LegionRuntime::Accessor::AccessorType::Generic, float> GfloatRA;
typedef LegionRuntime::Accessor::RegionAccessor < LegionRuntime::Accessor::AccessorType::Generic, RGBColor> GRGBColorRA;
typedef LegionRuntime::Accessor::RegionAccessor < LegionRuntime::Accessor::AccessorType::Generic, unsigned char> GucharRA;

//-----------------------------------------------------------------------------------------------------



#endif
