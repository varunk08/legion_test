/*
Ray tracer written using Legion
*/
#include <math.h>
#include <cmath>
#include <iostream>
#include "common_types.h"
#include "volumedata.h"
#include "xmlload.h"
#include "scene.h"
#include "boxobject.h" // remove later
#include "volume_geom.h"
#include "legion.h"

//-----------------------------------------------------------------------------------------------------
#define _USE_MATH_DEFINES;
using namespace std;
using namespace LegionRuntime::HighLevel;
using namespace LegionRuntime::Accessor;

//-----------------------------------------------------------------------------------------------------

//Function signatures
void GenerateRay(int x, int y, const Camera &camera, Ray &r, int xdim, int ydim);
LogicalRegion load_tf_data( int tf_size, cyColor* color_tf, float* alpha_tf, Context ctx, HighLevelRuntime *runtime);
LogicalRegion load_volume_data(int data_xdim, int data_ydim, int data_zdim, 
			       unsigned char *volume_data,
			       Context ctx, HighLevelRuntime *runtime);
bool init_volume_data ( unsigned char **volume_data, cyColor** color_tf, float** alpha_tf,
			unsigned char &minData, unsigned char &maxData, int &tf_size,
			int data_xdim, int data_ydim, int data_zdim,
			const char* data_file, const char* tf_file);
bool TraceSingleNode(const Ray &r, HitInfo &hInfo, const Node &node);
//-----------------------------------------------------------------------------------------------------


//Tasks
//Sub-tasks launched from top level task
void per_pixel_task ( const Task* task,
		     const std::vector<PhysicalRegion> &regions,
		     Context ctx, HighLevelRuntime *runtime)
{


  Ray* rays = (Ray*)task->local_args;
  //initialize node, object, transforms for each task
  Node *node = new Node;
  Sphere sphereObj;
  VolumeGeometry vol_obj(256, 256, 256, regions, ctx, runtime); //fix later - must get dimensions from args or regions
  PointLight *light = new PointLight();
  AmbientLight *amb_light = new AmbientLight();
  LightList lightList;
  MtlBlinn *mtlBlinn =  new MtlBlinn();

  //material init
  mtlBlinn->SetDiffuse(cyColor(1.0,0.5,0.5));
  mtlBlinn->SetSpecular(cyColor(0.7,0.7,0.7));
  //light init
  amb_light->SetIntensity(cyColor(0.2,0.2,0.2));
  light->SetIntensity(cyColor(0.7,0.7,0.7));
  light->SetPosition(Point3(0,10,0));
  lightList.push_back(amb_light);
  lightList.push_back(light);
  //node init
  node->SetName("only_node");
  node->SetObject(&vol_obj);
  //node->SetObject(&sphereObj);
  node->Scale(2,2,2);
  node->Rotate(Point3(0,0,1), 180);
  node->Rotate(Point3(0,1,0), 120);
  node->Translate (Point3(0,0,-5));
  node->SetObjTransform();
  node->SetMaterial(mtlBlinn);

  //init volume obj accessors
  GPoint3fRA acc_corner_pos = regions[1].get_field_accessor(FID_CORNER_POS).typeify<Point3f>();
  GucharRA acc_vol_data = regions[1].get_field_accessor(FID_VOL_DATA).typeify<unsigned char>();
  GPoint3fRA acc_gradients = regions[1].get_field_accessor(FID_GRADIENTS).typeify<Point3f>();
  GRGBColorRA acc_color_tf = regions[2].get_field_accessor(FID_COLOR_TF).typeify<RGBColor>();
  GfloatRA acc_alpha_tf = regions[2].get_field_accessor(FID_ALPHA_TF).typeify<float>();
  //1 - volume data; 2 - tf data
  vol_obj.init_logical_regions(task->regions[1].region, task->regions[2].region,
				acc_vol_data, acc_corner_pos, acc_gradients, acc_color_tf, acc_alpha_tf);
  vol_obj.init_tf_bounds(0, 255); //fix later - magic number
  vol_obj.set_lights(lightList);
  const int point = task->index_point.point_data[0];

  cout<<"Starting Rendering block: "<<point<<endl;
 
 //Writing data to logical region
  RegionAccessor<AccessorType::Generic, RGBColor> acc = regions[0].get_field_accessor(FID_OUT).typeify<RGBColor>();
  RGBColor col;
  col.r = 10.0f; col.g =10.0f; col.b=10.0f;
  Domain dom = runtime->get_index_space_domain(ctx, task->regions[0].region.get_index_space());
  Rect<1> rect = dom.get_rect<1>();
  int index = 0;

  for(GenericPointInRectIterator<1> pir(rect); pir; pir++)
    {
      Ray r = rays[index];
      HitInfo hInfo;
      hInfo.Init();

      //test for intersections
      if ( TraceSingleNode ( r, hInfo, *node)){

	//shade - fix magic number
	cyColor shade = hInfo.shade;
	float amp = 110.0f; //amp to control intensity - fix later
	col.r = amp * shade.r;
	col.g = amp * shade.g;
	col.b = amp * shade.b ;

      }

      //write data to logical region
      acc.write(DomainPoint::from_point<1>(pir.p), (RGBColor)col);
      col.r = 10.0f; col.g =10.0f; col.b=10.0f;
      ++index;
    }
  cout<<"Done rendering block "<<point<<endl;

}
//-----------------------------------------------------------------------------------------------------

//Volume render sub task called from per_pixel_task
/*
task arguments: shading parameters
logical regions: volume data
returns color as future
*/
void volume_render_task ( const Task* task,
			  const std::vector<PhysicalRegion> &regions,
			  Context ctx, HighLevelRuntime *runtime)
{
}

//-----------------------------------------------------------------------------------------------------
//Main task
void top_level_task( const Task* task,
		     const std::vector<PhysicalRegion> &regions,
		     Context ctx, HighLevelRuntime *runtime)
{


  //init objects used in scene
  Node         rootNode;
  Camera       camera;
  RenderImage  renderImage;
  MaterialList materials;
  LightList    lights;
  ObjFileList  objList;
  Sphere       theSphere;
  Plane        thePlane;
  BoxObject    theBoxObject;

  const char* filename = "scene.xml";
  /*write better version - most of the objects here are not needed except camera*/
    LoadScene(filename, rootNode, camera, renderImage,
	    materials, lights,  objList, theSphere,
	    theBoxObject,  thePlane );

  const char* data_file = "foot_8bit_256x256x256.raw";
  const char* tf_file = "foot_2.1dt";
  int data_xdim = 256; int data_ydim = 256; int data_zdim = 256;

  //Have to deallocate uc_vol_data, color_tf, alpha_tf later
  //  int size = data_xdim * data_ydim * data_zdim;
  unsigned char* uc_vol_data = NULL;
  cyColor* color_tf = NULL;
  float* alpha_tf = NULL;
  unsigned char minData; 
  unsigned char maxData;
  int tf_size;
  if( !init_volume_data (&uc_vol_data, &color_tf, &alpha_tf,
			 minData, maxData, tf_size,
			 data_xdim, data_ydim, data_zdim,
			 data_file, tf_file )){
    cout<<" Fatal error"<<endl;
    return;
  }

  LogicalRegion volume_lr = load_volume_data(data_xdim, data_ydim, data_zdim,
					     uc_vol_data, ctx, runtime);
  LogicalRegion tf_lr = load_tf_data(tf_size, color_tf, alpha_tf, ctx, runtime);
  
  delete [] uc_vol_data;
  delete [] color_tf;
  delete [] alpha_tf;

  //create rect of screen dimensions
  int xdim = 800;
  int ydim = 600;
  int size = xdim * ydim;
  Rect<1> elem_rect(Point<1>(0), Point<1>(size -1));

  //create index space
  IndexSpace is = runtime->create_index_space(ctx, Domain::from_rect<1>(elem_rect));

  //create output field space
  FieldSpace output_fs = runtime->create_field_space(ctx);
  {
    //field allocator for output field
    FieldAllocator allocator = runtime->create_field_allocator(ctx, output_fs);
    allocator.allocate_field(sizeof(RGBColor), FID_OUT);
  }

  //output logical region
  LogicalRegion output_lr = runtime->create_logical_region(ctx, is, output_fs);

  //Launch domain - our 2D grid
  Domain launch_domain = Domain::from_rect<1>(elem_rect);


  //Create Logical Partitions
  //!Make sure size is as multiple of num_regions!!
  int num_regions = 4;
  Rect<1> color_bounds(Point<1>(0), Point<1>(num_regions - 1));
  Domain color_domain = Domain::from_rect<1>(color_bounds);
  IndexPartition ip;
  int num_per_block = size / num_regions;
  Blockify<1> coloring ( num_per_block );
  ip = runtime->create_index_partition(ctx, is, coloring);
  LogicalPartition output_lp = runtime->get_logical_partition(ctx, output_lr, ip);

  cout<<"Creating Argument Map"<<endl;
  //arguments for each point in the 2D grid
  ArgumentMap arg_map;
  std::vector<Ray> rays;
  rays.clear();

  for(int i=0; i < size; ++i) {
    Ray r;
    int x = i % xdim;
    int y = i / xdim;
    GenerateRay(x, y, camera, r, xdim, ydim);
    rays.push_back(r);
  }

  for(int i=0; i< num_regions; ++i) {
    arg_map.set_point(DomainPoint::from_point<1>(Point<1>(i)),
		      TaskArgument(&rays[i * num_per_block], num_per_block * sizeof(Ray)));
  }

  cout<<"Launching tasks"<<endl;

  //Index launcher - for init values into logical region
  launch_domain = color_domain;
  IndexLauncher index_launcher ( PER_PIXEL_TASK_ID, 
				 launch_domain, 
				 TaskArgument(NULL,0), 
				 arg_map );
  //this example works without region requirements
  index_launcher.add_region_requirement( RegionRequirement( output_lp, 0,
							    WRITE_DISCARD,
							    EXCLUSIVE,
							    output_lr)); 
  index_launcher.add_region_requirement(RegionRequirement(volume_lr, 0, READ_ONLY, RELAXED, volume_lr));
  index_launcher.add_region_requirement(RegionRequirement(tf_lr, 0, READ_ONLY, RELAXED, tf_lr));

  //output 
  index_launcher.region_requirements[0].add_field(FID_OUT);

  //volume data
  index_launcher.region_requirements[1].add_field(FID_VOL_DATA);
  index_launcher.region_requirements[1].add_field(FID_CORNER_POS);
  index_launcher.region_requirements[1].add_field(FID_GRADIENTS);

  //transfer function data
  index_launcher.region_requirements[2].add_field(FID_COLOR_TF);
  index_launcher.region_requirements[2].add_field(FID_ALPHA_TF);

  runtime->execute_index_space(ctx, index_launcher);

  //output buffer
  RenderImage r_image;
  r_image.Init( xdim, ydim);

  //Launcher for reading physical region
  InlineLauncher output_launcher (RegionRequirement (output_lr, READ_ONLY, EXCLUSIVE, output_lr));
  output_launcher.requirement.add_field(FID_OUT);
  PhysicalRegion output_pr = runtime->map_region(ctx, output_launcher);
  output_pr.wait_until_valid();
  RegionAccessor<AccessorType::Generic, RGBColor> acc = output_pr.get_field_accessor(FID_OUT).typeify<RGBColor>();
  int index = 0;

  //read from logical region and write to buffer renderimage
  for (GenericPointInRectIterator<1> pir (elem_rect); pir; pir++) {
    RGBColor val = (RGBColor)acc.read(DomainPoint::from_point<1>(pir.p));
    Color24 col;
    col.r  = val.r;
    col.g  = val.g;
    col. b = val.b;
    r_image.PutPixel(index, col, val.r);
    index++;
  }
  r_image.SaveImage("renderimage.ppm");

  //Destroy all - index spaces, field spaces, logical regions
  runtime->destroy_logical_region( ctx, output_lr);
  runtime->destroy_field_space( ctx, output_fs);
  runtime->destroy_index_space( ctx, is);
}
//-----------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{

  HighLevelRuntime::set_top_level_task_id( TOP_LEVEL_ID);

  HighLevelRuntime::register_legion_task<top_level_task>(TOP_LEVEL_ID,
							 Processor::LOC_PROC,
							 true,
							 false );
  HighLevelRuntime::register_legion_task<per_pixel_task>(PER_PIXEL_TASK_ID,
							 Processor::LOC_PROC,
							 true,
							 false);
  HighLevelRuntime::register_legion_task<volume_render_task>(VOLUME_RENDER_TASK_ID,
							     Processor::LOC_PROC,
							     true,
							     false);

  return HighLevelRuntime::start(argc, argv);

}
//-----------------------------------------------------------------------------------------------------
bool init_volume_data ( unsigned char **volData, cyColor **color_tf, float **alpha_tf,
			unsigned char &minData, unsigned char &maxData, int &tf_size,
			int data_xdim, int data_ydim, int data_zdim,
			const char* data_file, const char* tf_file)
{

  VolumeData DataLoader;

  if( DataLoader.Load(data_file, data_xdim, data_ydim, data_zdim, volData) 
      && DataLoader.LoadTF(tf_file) ) {
    std::cout<<"Loaded data file: "<<data_file<<std::endl;
    DataLoader.GetTransferFunction( color_tf, alpha_tf, tf_size, minData, maxData );
    //theBoxObject.SetTransferFunction(color_tf, alpha_tf, tf_size, minData, maxData);
  }
  else {
    cout<<"Failed loading volume data file"<<endl;
    return false;
  }
  return true;
}
int getIndex(int x, int y, int z, int xdim, int ydim) 
{
  return ( z * ydim + y ) * xdim + x;
}
LogicalRegion load_tf_data( int tf_size, cyColor* color_tf, float* alpha_tf, Context ctx, HighLevelRuntime *runtime)
{
  //Index space and field space
  Rect<1> tf_pts_rect(Point<1>(0), Point<1>(tf_size-1));
  IndexSpace tf_is = runtime->create_index_space( ctx, Domain::from_rect<1>(tf_pts_rect));
  FieldSpace tf_fs = runtime->create_field_space(ctx);

  FieldAllocator tf_falloc = runtime->create_field_allocator(ctx, tf_fs);
  tf_falloc.allocate_field ( sizeof ( Point3f ) , FID_COLOR_TF);
  tf_falloc.allocate_field ( sizeof (float), FID_ALPHA_TF);

  //Logical region
  LogicalRegion tf_lr = runtime->create_logical_region( ctx, tf_is, tf_fs);
  runtime->attach_name( tf_lr, "transfer_lr");

  //Region requirements, accessor and writing data for transfer functions
  RegionRequirement tf_req(tf_lr, READ_WRITE, EXCLUSIVE, tf_lr);
  tf_req.add_field(FID_COLOR_TF);
  tf_req.add_field(FID_ALPHA_TF);

  PhysicalRegion tf_pr = runtime->map_region(ctx, tf_req);
  tf_pr.wait_until_valid();

  RegionAccessor< AccessorType::Generic, RGBColor> fa_tf_color = tf_pr.get_field_accessor(FID_COLOR_TF).typeify<RGBColor>();
  RegionAccessor<AccessorType::Generic, float> fa_tf_alpha = tf_pr.get_field_accessor(FID_ALPHA_TF).typeify<float>();

  //write data in color_tf and alpha_tf to logical region
  cout<<"Writing transfer function to logical region"<<endl;
  GenericPointInRectIterator<1> tf_itr(tf_pts_rect);
  for(int i=0; i<tf_size; ++i, tf_itr++){
    RGBColor col;
    col.r = color_tf[i].r; col.g = color_tf[i].g; col.b = color_tf[i].b;
    fa_tf_color.write(DomainPoint::from_point<1>(tf_itr.p), col);
    fa_tf_alpha.write(DomainPoint::from_point<1>(tf_itr.p), alpha_tf[i]);
  }

  runtime->unmap_region(ctx, tf_pr);

  return tf_lr;

}
LogicalRegion load_volume_data(int data_xdim, int data_ydim, int data_zdim,
			       unsigned char *volumeData,
			       Context ctx, HighLevelRuntime *runtime)
{

  int num_points = data_xdim * data_ydim * data_zdim;

  //Index space
  Rect<1> vol_points_rect(Point<1>(0), Point<1>(num_points-1));
  IndexSpace volume_data_is = runtime->create_index_space( ctx,Domain::from_rect<1>(vol_points_rect));

  //Field space
  FieldSpace volume_data_fs = runtime->create_field_space(ctx);


  FieldAllocator falloc = runtime->create_field_allocator(ctx, volume_data_fs);
  falloc.allocate_field(sizeof(unsigned char), FID_VOL_DATA);
  falloc.allocate_field(sizeof(Point3f), FID_CORNER_POS);
  falloc.allocate_field(sizeof(Point3f), FID_GRADIENTS);

  //Make logical region
  LogicalRegion volume_data_lr = runtime->create_logical_region(ctx, volume_data_is, volume_data_fs);
  runtime->attach_name(volume_data_lr, "volume_data_lr");

  //Create Region requirement
  RegionRequirement vol_req (volume_data_lr, READ_ONLY,RELAXED , volume_data_lr);
  vol_req.add_field(FID_VOL_DATA);
  vol_req.add_field(FID_CORNER_POS);
  vol_req.add_field(FID_GRADIENTS);

  //Create Physical regions
  PhysicalRegion vol_data_pr = runtime->map_region(ctx, vol_req);

  //Create Region accessors
  vol_data_pr.wait_until_valid();
  RegionAccessor <AccessorType::Generic, unsigned char> fa_vol_data = vol_data_pr.get_field_accessor(FID_VOL_DATA).typeify<unsigned char>();
  RegionAccessor<AccessorType::Generic, Point3f> fa_corner_pos = vol_data_pr.get_field_accessor(FID_CORNER_POS).typeify<Point3f>();
  RegionAccessor<AccessorType::Generic, Point3f> fa_gradients = vol_data_pr.get_field_accessor(FID_GRADIENTS).typeify<Point3f>();

  //Init logical region with data from VolumeData to (1) FID_VOL_DATA (2) FID_CORNER_POS (3) FID_GRADIENTS
  //write volume data to logical region
   GenericPointInRectIterator<1> itr(vol_points_rect);
  for(int i=0; i<num_points; ++i, itr++){
    fa_vol_data.write(DomainPoint::from_point<1>(itr.p), volumeData[i]);
  }

  //Calculate corner positions
  cout<<"Writing volume data to logical region"<<endl;
  float startX = -1.0f;  float startY = -1.0f;  float startZ = -1.0f;
  float newX = 0.0f, newY = 0.0f, newZ = 0.0f;
  Point3f c_pos;
  GenericPointInRectIterator<1> cpos_itr(vol_points_rect);
  for(int z = 0; z < data_zdim; z++)
    {
      newZ = startZ + z * (2.0 / (data_zdim - 1));
      for(int y = 0; y < data_ydim; y++)
	{
	  newY = startY + y * (2.0 / (data_ydim - 1));
	  for(int x = 0; x < data_xdim; x++)
	    {
	      newX = startX + x * (2.0 / (data_xdim - 1));
	      //Write corner position info to region

	      c_pos.x = newX; c_pos.y = newY; c_pos.z = newZ;

	      //write corner position info to logical region
	      fa_corner_pos.write(DomainPoint::from_point<1>(Point<1>(getIndex(x,y,z,data_xdim,data_ydim))),
				  (Point3f) c_pos);
	      
	      //calculate gradients for each point
	      float xn, xp, yn, yp, zn, zp;
	      if(x > 0) xp = (uint)volumeData[getIndex(x-1,y,z,data_xdim, data_ydim)];
	      else xp=0;
	      if(x < data_xdim-1) xn = (uint)volumeData[getIndex(x+1,y,z,data_xdim, data_ydim)];
	      else xn = 0;
	      if(y > 0) yp = (uint)volumeData[getIndex(x,y-1,z,data_xdim, data_ydim)];
	      else  yp = 0;
	      if(y < data_ydim-1) yn = (uint)volumeData[getIndex(x,y+1,z,data_xdim, data_ydim)];
	      else yn = 0;
	      if(z > 0) zp = (uint)volumeData[getIndex(x,y,z-1,data_xdim, data_ydim)];
	      else zp=0;
	      if(z < data_zdim-1) zn = (uint)volumeData[getIndex(x,y,z+1,data_xdim, data_ydim)];
	      else zn = 0;
	      Point3f N; 
	      N.x = xn - xp;
	      N.y =  yn - yp;
	      N.z = zn - zp;

	      //write gradients to logical region
	      fa_gradients.write(DomainPoint::from_point<1>(Point<1>(getIndex(x,y,z,data_xdim,data_ydim))), N);
	      cpos_itr++;

	    }
	}
    }
  runtime->unmap_region(ctx, vol_data_pr);

  return volume_data_lr;
}
//-----------------------------------------------------------------------------------------------------
void GenerateRay(int x, int y, const Camera &camera, Ray &r, int xdim, int ydim)
{
  float alpha = camera.fov;
  float l = 1.0f;
  float h = l * tan (alpha / 2.0 * (M_PI/180));
  float aspectRatio = (float) xdim / ydim;
  float w = aspectRatio * abs(h);
  float dx = (2 * abs(w)) / xdim;
  float dy = -(2 * abs(h)) / ydim;
  float dxx = dx/2, dyy = dy/2;
  Point3 K(-w, h, -l);


  K.x += x * dx;
  K.y += y * dy;
  K.x += dxx;
  K.y += dyy;

  Matrix3 RotMat;
  cyPoint3f f = camera.dir;
  f. Normalize();

  cyPoint3f s = f.Cross(camera.up);
  s. Normalize();

  cyPoint3f u = s.Cross(f);
  u.Normalize();
  const float pts[9] = {s.x, u.x, -f.x,
			s.y, u.y, -f.y,
			s.z, u.z, -f.z  };
  RotMat.Set(pts);
  Point3 cam_pos(camera.pos);
  r.p = cam_pos;
  r.dir =  K * RotMat;

  r.dir.Normalize();
}
//-----------------------------------------------------------------------------------------------------
bool Box::IntersectRay(const Ray &ray, float t_max) const
{
    Ray r = ray;
    float tmin  =  -t_max;
    float tmax = t_max;
    //if ray is inside box - return true
    if (IsInside(r.p)) return true;
    //get pairs of planes - x , y, z
    // 0:(x_min,y_min,z_min), 1:(x_max,y_min,z_min)
    // 2:(x_min,y_max,z_min), 3:(x_max,y_max,z_min)
    // 4:(x_min,y_min,z_max), 5:(x_max,y_min,z_max)
    // 6:(x_min,y_max,z_max), 7:(x_max,y_max,z_max)  }  }
    float xl = Corner(0).x;
    float xh = Corner(3).x;
    float yl = Corner(0).y;
    float yh = Corner(2).y;
    float zl = Corner(0).z;
    float zh = Corner(5).z;
    
    //Check intersection for X planes
    if(r.p.x != 0.0){
        float tx1 = (xl - ray.p.x)/ ray.dir.x;
        float tx2 = (xh - ray.p.x)/ ray.dir.x;
        
        tmin = utils::max((float)tmin,(float) utils::min((float)tx1,tx2));
        tmax = utils::min((float)tmax,(float) utils::max((float)tx1, tx2));
    }
    //Check intersection for Y planes
    if(r.p.y != 0.0){
        float tx1 = (yl - ray.p.y)/ ray.dir.y;
        float tx2 = (yh - ray.p.y)/ ray.dir.y;
        
        tmin = max(tmin, min(tx1,tx2));
        tmax = min(tmax, max(tx1, tx2));
    }
    
    //Check intersection for Z planes
    if(r.p.z != 0.0){
        float tx1 = (zl - ray.p.z)/ ray.dir.z;
        float tx2 = (zh - ray.p.z)/ ray.dir.z;
        
        tmin = max(tmin, min(tx1,tx2));
        tmax = min(tmax, max(tx1, tx2));
    }
    
    return tmax>=tmin;
    
}
//-----------------------------------------------------------------------------------------------------
bool TraceSingleNode(const Ray &r, HitInfo &hInfo, const Node &node)
{
    Ray ray = r;
    ray = node.ToNodeCoords(r);

    const Object* obj = node.GetObject();
    bool objHitTest = false;

    
    if(obj)    {
        if(obj->IntersectRay(ray, hInfo)){
            objHitTest=true;
            hInfo.node = &node;
        }
    }
    
    if(objHitTest){
       node.FromNodeCoords(hInfo);
    }
    return objHitTest;
}
//-----------------------------------------------------------------------------------------------------
