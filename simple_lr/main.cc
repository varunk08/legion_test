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
}
void top_level_task( const HighLevel::Task* task,
		     const std::vector<HighLevel::PhysicalRegion> &regions,
		     HighLevel::Context ctx, HighLevel::HighLevelRuntime *runtime)
{

  int xdim = 400;
  int ydim = 400;
  cout<< "Entering top level task"<<endl;
  HighLevel::Rect<1> 

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
