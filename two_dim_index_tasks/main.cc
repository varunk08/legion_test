/*
Program to create logical region for a grid of points and perform some process on them

 */
#include <iostream>
#include "legion.h"

using namespace LegionRuntime;
using namespace std;

enum TaskID{
  TOP_LEVEL_ID,
  PER_POINT_TASK_ID,
};

int per_point_task ( const HighLevel::Task* task,
		     const std::vector<HighLevel::PhysicalRegion> &regions,
		     HighLevel::Context ctx, HighLevel::HighLevelRuntime *runtime)
{
  int input = *((const int*) task->local_args);
  return (1 + input);
}
void top_level_task( const HighLevel::Task* task,
		     const std::vector<HighLevel::PhysicalRegion> &regions,
		     HighLevel::Context ctx, HighLevel::HighLevelRuntime *runtime)
{

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

  HighLevel::Domain launch_domain = HighLevel::Domain::from_rect<2>(elem_rect);
  HighLevel::ArgumentMap arg_map;
  for(int i = 0; i < xdim; ++i) {
    for(int j=0; j < ydim; ++j) {
      int input = i * j;
      int pt[2];    pt[0] =i; pt[1] = j;
      arg_map.set_point( HighLevel::DomainPoint::from_point<2>(Point<2>(pt)) , HighLevel::TaskArgument (&input, sizeof(input)));
    }
  }

  HighLevel::IndexLauncher index_launcher ( PER_POINT_TASK_ID, launch_domain, HighLevel::TaskArgument(NULL,0), arg_map);
  HighLevel::FutureMap fm = runtime->execute_index_space(ctx, index_launcher);
  fm.wait_all_results();
  for(int i = 0; i < xdim; ++i) {
    for(int j =0; j< ydim; ++j) {
      int pt[2];  pt[0] =i; pt[1] = j;

      int received = fm.get_result<int>( HighLevel::DomainPoint::from_point<2>(Point<2>(pt)));
      cout<<"("<<i<<","<<j<<"): "<<received<<" ";

    }
    cout<<"\n";
  }
}

int main(int argc, char **argv)
{
  HighLevel::HighLevelRuntime::set_top_level_task_id( TOP_LEVEL_ID);

  HighLevel::HighLevelRuntime::register_legion_task<top_level_task>(TOP_LEVEL_ID,
								    HighLevel::Processor::LOC_PROC,
								    true,
								    false );
  HighLevel::HighLevelRuntime::register_legion_task<int, per_point_task>(PER_POINT_TASK_ID,
									 HighLevel::Processor::LOC_PROC,
									 true,
									 false);
  return HighLevel::HighLevelRuntime::start(argc, argv);

}
