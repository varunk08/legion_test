/*
Ray tracer written using Legion
*/
#include <math.h>
#include <cmath>
#include <iostream>
#include "legion.h"
#include "xmlload.h"
#include "scene.h"
#include "render_object.h"
#define _USE_MATH_DEFINES;
using namespace LegionRuntime::HighLevel;
using namespace LegionRuntime::Accessor;
using namespace std;

enum TaskID{
  TOP_LEVEL_ID,
  PER_POINT_TASK_ID,
  CHECK_TASK_ID,
};
enum FieldSpaceID{
  FID_IN,
  FID_OUT,
};

//Info needed by each task to perform rendering calculations
struct SceneData 
{
  Ray ray; //ray info
  int coords[2];
};

void GenerateRay(int x, int y, const Camera &camera, Ray &r)
{
  float alpha = camera.fov;
  float l = 1.0f;
  float h = l * tan (alpha / 2.0 * (M_PI/180));
  float aspectRatio = (float) camera.imgWidth / camera.imgHeight;
  float w = aspectRatio * abs(h);
  float dx = (2 * abs(w)) / camera.imgWidth;
  float dy = -(2 * abs(h)) / camera.imgHeight;
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
  const float pts[9] = {s.x, u.x , -f.x,
		       s.y, u.y, -f.y,
		       s.z, u.z, -f.z  };
  RotMat.Set(pts);
  Point3 cam_pos(camera.pos);
  r.p = cam_pos;
  r.dir =  K;
  r.dir =  r.dir * RotMat;
  r.dir.Normalize();



}
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
void check_task ( const Task* task,
		  const std::vector<PhysicalRegion> &regions,
		  Context ctx, HighLevelRuntime *runtime)
{

  RegionAccessor<AccessorType::Generic, float> acc = regions[0].get_field_accessor(FID_IN).typeify<float>();
  Domain dom = runtime->get_index_space_domain(ctx,
					       task->regions[0].region.get_index_space());
  Rect<2> rect = dom.get_rect<2>();
  for (GenericPointInRectIterator<2> pir (rect); pir; pir++) {
    float val = acc.read(DomainPoint::from_point<2>(pir.p));
    cout<<"val: "<<val<<endl;
  }

}

bool TraceSingleNode(const Ray &r, HitInfo &hInfo, const Node &node)
{

    Ray ray = r;
    ray = node.ToNodeCoords(r);

    const Object* obj = node.GetObject();
    bool objHitTest = false;
    bool childHit = false;
    
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

//Sub-tasks launched from top level task
void per_point_task ( const Task* task,
		     const std::vector<PhysicalRegion> &regions,
		     Context ctx, HighLevelRuntime *runtime)
{


  /*
1. init node
2. init materials 
3. init transformations
4. perform ray tracing
5. write to logical region
*/

  /*SceneData sd  = *((SceneData*)task->local_args);

  //initialize node, object, transforms for each task
  Node* node = new Node;
  Sphere sphereObj;
  node->SetObject(&sphereObj);
  node->Scale(2,2,2);
  node->Rotate(Point3(0,1,0), 30);
  node->Translate (Point3(0,0,-1));

  Ray r = sd.ray;
  HitInfo hInfo;
  hInfo.Init();
  float hit = 0.0f;
  if ( TraceSingleNode ( r, hInfo, *node)){
    hit = 1.0f;
  }
  */
  FieldID fid = *(task->regions[0].privilege_fields.begin());
  const int point = task->index_point.point_data[0];
  
  RegionAccessor<AccessorType::Generic, float> acc = 
    regions[0].get_field_accessor(fid).typeify<float>();

  Domain dom = runtime->get_index_space_domain(ctx,
					       task->regions[0].region.get_index_space());
  Rect<1> rect = dom.get_rect<1>();
  for(GenericPointInRectIterator<1> pir(rect); pir; pir++)
    {
      acc.write(DomainPoint::from_point<1>(pir.p), 0.5f);
    }


  /*//TEST */
    int* input = ((int*)task->local_args);
    cout<<"Test: "<<input[0]<<" "<<input[1]<<endl;


}

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

  LoadScene(filename, rootNode, camera, renderImage,
	    materials, lights,  objList, theSphere,
	    theBoxObject,  thePlane );

  //init render object
  /*  RenderObject render_object;
  Node* child_node = rootNode.GetChild(0);
  render_object.SetObject( *((Sphere *)child_node->GetObject()) ); //get the first child
  render_object.SetMaterial( child_node->GetMaterial() );
  //render_object.SetLightList(lights);
  render_object.SetTransformation(child_node->GetTransform());
  */

  //create rect of screen dimensions
  int xdim = 800;
  int ydim = 600;
  int size = xdim * ydim;
  //  int lo[2] = {0,0}; 
  //  int hi[2] = {xdim -1, ydim -1};
  //  Point<2> low(lo);  Point<2> high(hi);
  //  Rect<2> elem_rect(low,high); 
  Rect<1> elem_rect(Point<1>(0), Point<1>(size -1));

  //create index space
  IndexSpace is = runtime->create_index_space(ctx, Domain::from_rect<1>(elem_rect));

  //create input field space
  FieldSpace input_fs = runtime->create_field_space(ctx);
  {
    //field space allocator
    FieldAllocator allocator = runtime->create_field_allocator(ctx, input_fs);
    allocator.allocate_field(sizeof(float), FID_IN);
  }
  //create output field space
  FieldSpace output_fs = runtime->create_field_space(ctx);
  {
    //field allocator for output field
    FieldAllocator allocator = runtime->create_field_allocator(ctx, output_fs);
    allocator.allocate_field(sizeof(float), FID_OUT);
  }

  //create logical regions for input and output
  LogicalRegion input_lr = runtime->create_logical_region(ctx, is, input_fs);

  //output logical region
  LogicalRegion output_lr = runtime->create_logical_region(ctx, is, output_fs);

  //Launch domain - our 2D grid
  Domain launch_domain = Domain::from_rect<1>(elem_rect);


  //Create Logical Partitions
  int num_regions = 48;
  Rect<1> color_bounds(Point<1>(0), Point<1>(num_regions - 1));
  Domain color_domain = Domain::from_rect<1>(color_bounds);
  IndexPartition ip;
  Blockify<1> coloring ( size / num_regions );
  ip = runtime->create_index_partition(ctx, is, coloring);
  LogicalPartition output_lp = runtime->get_logical_partition(ctx, output_lr, ip);


  //arguments for each point in the 2D grid
  ArgumentMap arg_map;

  for(int i=0; i < size; ++i) {
    int x = i % xdim;
    int y = i / xdim;
    int pt[2] = { x, y };
    //    Ray prim_ray;
    //    GenerateRay(x, y, camera, prim_ray);
    //    SceneData scene_data;
    //    scene_data.ray = prim_ray;
    //    scene_data.rootNode = rootNode;
    //    scene_data.output_lr = output_lr;
    //    scene_data.coords[0] = x;      scene_data.coords[1] = y;
    //arg_map.set_point( DomainPoint::from_point<2>(Point<2>(pt)), TaskArgument (&scene_data, sizeof(SceneData)));
    arg_map.set_point( DomainPoint::from_point<1>(Point<1>(i)), TaskArgument (&pt, 2 * sizeof(int)));
  }

  cout<<"Launching tasks"<<endl;
  //Index launcher - for init values into logical region
  IndexLauncher index_launcher ( PER_POINT_TASK_ID, 
				 color_domain, 
				 TaskArgument(NULL,0), 
				 arg_map );
  //this example works without region requirements
  index_launcher.add_region_requirement( RegionRequirement( output_lp, 0,
							    WRITE_DISCARD,
							    EXCLUSIVE,
							    output_lr));
  index_launcher.region_requirements[0].add_field(FID_OUT);
  runtime->execute_index_space(ctx, index_launcher);

  float zbuffer[xdim * ydim];
  RenderImage r_image;
  r_image.Init( xdim, ydim);

  InlineLauncher output_launcher (RegionRequirement (output_lr, READ_ONLY, EXCLUSIVE, output_lr));
  output_launcher.requirement.add_field(FID_OUT);
  PhysicalRegion output_pr = runtime->map_region(ctx, output_launcher);

  RegionAccessor<AccessorType::Generic, float> acc = output_pr.get_field_accessor(FID_OUT).typeify<float>();
  int index = 0;

  for (GenericPointInRectIterator<1> pir (elem_rect); pir; pir++) {
    float val = acc.read(DomainPoint::from_point<1>(pir.p));
    zbuffer[index] = val;
    Color24 col;  col.r = val*150; col.g = val * 150; col.b = 0.0f;
    r_image.PutPixel(index,col, val);
    index++;
  }
  r_image.SaveImage("renderimage.ppm");


  //Task Launcher for check tasks
  /*  TaskLauncher check_launcher( CHECK_TASK_ID, TaskArgument(NULL,0));
  check_launcher.add_region_requirement( RegionRequirement(input_lr,
							    READ_ONLY,
							    EXCLUSIVE,
							    input_lr));
  check_launcher.region_requirements[0].add_field ( FID_IN );
  runtime->execute_task(ctx, check_launcher);
  */  
  //Destroy all - index spaces, field spaces, logical regions
  runtime->destroy_logical_region( ctx, input_lr );
  runtime->destroy_logical_region( ctx, output_lr);
  runtime->destroy_field_space( ctx, input_fs);
  runtime->destroy_field_space( ctx, output_fs);
  runtime->destroy_index_space( ctx, is);
}

int main(int argc, char **argv)
{
  HighLevelRuntime::set_top_level_task_id( TOP_LEVEL_ID);

  HighLevelRuntime::register_legion_task<top_level_task>(TOP_LEVEL_ID,
								    Processor::LOC_PROC,
								    true,
								    false );
  HighLevelRuntime::register_legion_task<per_point_task>(PER_POINT_TASK_ID,
									 Processor::LOC_PROC,
									 true,
									 false);
  HighLevelRuntime::register_legion_task<check_task>( CHECK_TASK_ID,
								 Processor::LOC_PROC,
								 true,
								 false);

  return HighLevelRuntime::start(argc, argv);

}
