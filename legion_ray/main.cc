/*
Ray tracer written using Legion
*/
#include <iostream>
#include "legion.h"
#include "xmlload.h"
#include "scene.h"

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
    // 6:(x_min,y_max,z_max), 7:(x_max,y_max,z_max)
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
float GenLight::Shadow(Ray ray, float t_max)
{
    return 1.0; //direct
}

//Sub-tasks launched from top level task
void per_point_task ( const Task* task,
		     const std::vector<PhysicalRegion> &regions,
		     Context ctx, HighLevelRuntime *runtime)
{
  const int* input = ((const int*) task->local_args);
  int i = input[0]; int j = input[1];
  //TODO - write argument to logical region
  FieldID fid = *(task->regions[0].privilege_fields.begin());
  const int point = task->index_point.point_data[0];
  RegionAccessor<AccessorType::Generic, double> acc = regions[0].get_field_accessor(fid).typeify<double>();
  Domain dom = runtime->get_index_space_domain(ctx,
					       task->regions[0].region.get_index_space());
  Rect<2> rect = dom.get_rect<2>();
  int index = (5 * j) + i;
  cout<<"writing: "<<i<<" , "<<j<<"--->"<<index<<endl;
  acc.write(DomainPoint::from_point<2>(Point<2>(input)), index);

}
void check_task ( const Task* task,
		  const std::vector<PhysicalRegion> &regions,
		  Context ctx, HighLevelRuntime *runtime)
{

  RegionAccessor<AccessorType::Generic, double> acc = regions[0].get_field_accessor(FID_IN).typeify<double>();
  Domain dom = runtime->get_index_space_domain(ctx,
					       task->regions[0].region.get_index_space());
  Rect<2> rect = dom.get_rect<2>();
  for (GenericPointInRectIterator<2> pir (rect); pir; pir++) {
    double val = acc.read(DomainPoint::from_point<2>(pir.p));
    cout<<"val: "<<val<<endl;
  }

}

//Main task
void top_level_task( const Task* task,
		     const std::vector<PhysicalRegion> &regions,
		     Context ctx, HighLevelRuntime *runtime)
{
  //load xml file
  /* Have all global objects here - create logical regions for them as required*/
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

  LoadScene(filename, 
	    rootNode,
	    camera,
	    renderImage,
	    materials,
	    lights,
	    objList,
	    theSphere,
	    theBoxObject,
	    thePlane );

  //init objects used in scene
  //for each pixel spawn task for recursive ray-tracing of the scene graph
  //return future map of colors or simple write color into logical regions - simple
  int xdim = 5;
  int ydim = 5;
  int lo[2]; 
  int hi[2];
  lo[0] = 0; lo[1] = 0;
  hi[0] = xdim -1; hi[1] = ydim-1;
  Point<2> low(lo);
  Point<2> high(hi);
  cout<< "Entering top level task"<<endl;
  Rect<2> elem_rect(low,high); //(0,0) to (xdim-1,ydim-1)

  //create index space
  IndexSpace is = runtime->create_index_space(ctx, Domain::from_rect<2>(elem_rect));
  runtime->attach_name(is, "is");

  //create input field space
  FieldSpace input_fs = runtime->create_field_space(ctx);
  runtime->attach_name(input_fs, "input_fs");
  //field space allocator
  FieldAllocator allocator = runtime->create_field_allocator(ctx, input_fs);
  allocator.allocate_field(sizeof(double), FID_IN);
  runtime->attach_name(input_fs, FID_IN, "in");

  //create output field space
  FieldSpace output_fs = runtime->create_field_space(ctx);
  runtime->attach_name( output_fs, "out");
  //field allocator for out put field
  allocator = runtime->create_field_allocator(ctx, output_fs);
  allocator.allocate_field(sizeof(double), FID_OUT);
  runtime->attach_name(output_fs, FID_OUT, "out");

  //create logical regions for input and output
  LogicalRegion input_lr = runtime->create_logical_region(ctx, is, input_fs);
  runtime->attach_name( input_lr, "input_lr");
  LogicalRegion output_lr = runtime->create_logical_region(ctx, is, output_fs);
  runtime->attach_name( output_lr, "output_lr");

  //Launch domain - our 2D grid
  Domain launch_domain = Domain::from_rect<2>(elem_rect);

  //arguments for each point in the 2D grid
  ArgumentMap arg_map;
  
  for(int i = 0; i < xdim; ++i) {
    for(int j=0; j < ydim; ++j) {
      int pt[2];    pt[0] =i; pt[1] = j;
      arg_map.set_point( DomainPoint::from_point<2>(Point<2>(pt)) , 
			 TaskArgument (&pt, 2 * sizeof(int)));
    }
  }
  
  //Index launcher - for init values into logical region
  IndexLauncher index_launcher ( PER_POINT_TASK_ID, 
					    launch_domain, 
					    TaskArgument(NULL,0), 
					    arg_map );

  index_launcher.add_region_requirement( RegionRequirement( input_lr, 
							    WRITE_DISCARD,
							    EXCLUSIVE,
							    input_lr));
  index_launcher.region_requirements[0].add_field(FID_IN);
  runtime->execute_index_space(ctx, index_launcher);
  /*					
  HighLevel::FutureMap fm = runtime->execute_index_space(ctx, index_launcher);
  fm.wait_all_results();
  for(int i = 0; i < xdim; ++i) {
    for(int j =0; j< ydim; ++j) {
      int pt[2];  pt[0] =i; pt[1] = j;
      int received = fm.get_result<int>( HighLevel::DomainPoint::from_point<2>(Point<2>(pt)));
      cout<<"("<<i<<","<<j<<"): "<<received<<" ";
    }
    cout<<"\n";
    }*/

  //Task Launcher for check tasks
   TaskLauncher check_launcher( CHECK_TASK_ID, TaskArgument(NULL,0));
  check_launcher.add_region_requirement( RegionRequirement(input_lr,
							    READ_ONLY,
							    EXCLUSIVE,
							    input_lr));
  check_launcher.region_requirements[0].add_field ( FID_IN );
  runtime->execute_task(ctx, check_launcher);
  
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
