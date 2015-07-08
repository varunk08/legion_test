/*
Program to create logical region for a grid of points and perform some process on them

 */
#include <iostream>
#include "legion.h"

using namespace LegionRuntime::HighLevel;
using namespace LegionRuntime::Accessor;
using namespace std;

enum TaskID{
  TOP_LEVEL_ID,
  TEST_TASK_ID,
};
enum FieldSpaceID{
  FID_IN,
  FID_OUT,
};

struct TestData
{
  int a;
  float b;
  double c;
  char* name;
};
void test_task( const Task* task,
		     const std::vector<PhysicalRegion> &regions,
		     Context ctx, HighLevelRuntime *runtime)
{
  printf("\nExecuting child task\n");
  TestData input = *((TestData *)task->args);
  printf("a: %d, b: %f, c: %f, name: %s",input.a,input.b,input.c, input.name);
}
//Main task
void top_level_task( const Task* task,
		     const std::vector<PhysicalRegion> &regions,
		     Context ctx, HighLevelRuntime *runtime)
{

  //sample data
  TestData tda, tdb;
  tda.a = 1; tda.b = 10.0f; tda.c = 100.123f; tda.name = "task-1";
  tdb.a = 2; tdb.b = 20.0f; tdb.c = 200.123f; tdb.name = "task-2";


  TaskLauncher launcher(TEST_TASK_ID, TaskArgument(&tda, sizeof(TestData)));
  TaskLauncher launcher2(TEST_TASK_ID, TaskArgument(&tdb, sizeof(TestData)));
  runtime->execute_task(ctx, launcher);
  runtime->execute_task(ctx,launcher2);


}

int main(int argc, char **argv)
{
  HighLevelRuntime::set_top_level_task_id( TOP_LEVEL_ID);

  HighLevelRuntime::register_legion_task<top_level_task>(TOP_LEVEL_ID,
							 Processor::LOC_PROC,
							 true,
							 false );


  HighLevelRuntime::register_legion_task<test_task>(TEST_TASK_ID,
						    Processor::LOC_PROC,
						    true,
						    false );



  return HighLevelRuntime::start(argc, argv);

}
