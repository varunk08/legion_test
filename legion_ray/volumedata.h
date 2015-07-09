#ifndef VOLUMEDATA_H
#define VOLUMEDATA_H

/* Class for loading and parsing volume data; analyzing histogram and creating transfer functions */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "cyColor.h"
#include "cyPoint.h"
#include "utils.h"

class VolumeData
{

 public:
  int tf_size;
  int n_bins;
  int xdim, ydim, zdim;
  cyColor* colortf;
  float* alphatf;
  cyPoint4f* tf_rgba;
  unsigned char maxData, minData;
  //constructor
  VolumeData() {};

  /*Parses volume data*/
  void CreateHistogram(uchar* volume_data, int size);

  /*Loads volume data from file; calls CreateHistogram; calls CreateTransferFunction after parsing data*/
  bool Load(const char* filename, int xdim, int ydim, int zdim, uchar** data);

  /*Loads the transer function file; ImageVis3D format .1dt*/
  bool LoadTF(const char* filename);
  
  /*Creates transfer function;*/
  void CreateTransferFunction();

  /*Getter for transfer functions*/
  void GetTransferFunction(cyColor** color_tf, float** alpha_tf, int &tf_size,unsigned char &min, unsigned char &max);

};

#endif
