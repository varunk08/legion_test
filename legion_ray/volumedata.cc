#include "volumedata.h"
bool VolumeData::Load(const char* filename, int xdim, int ydim, int zdim, uchar** data)
{

  FILE *pFile;
  const int size = xdim*ydim*zdim;
  
  pFile = fopen(filename, "rb");
  if(NULL == pFile)
    {
      return false;
    }
  uchar *pVolume = new uchar[size];
  fread(pVolume,sizeof(unsigned char),size,pFile);
  fclose(pFile);
  CreateHistogram(pVolume, size);

  *data = pVolume;
  return true;

}

bool VolumeData::LoadTF(const char* filename)
{
   
  std::ifstream infile(filename, std::ifstream::in);
  std::string instr;
  if( !infile.good()){
    std::cout<<"Unable to load "<<filename<<std::endl;
    return false;
  }

  std::getline(infile, instr);
  const int num = atoi ( instr.c_str() );
  std::cout<<"Data: "<< num <<std::endl;
  tf_rgba = new cyPoint4f[num];
  this->colortf = new cyColor[num];
  this->alphatf = new float[num];

  int i = 0;

 
  while(std::getline(infile,instr)){
   
    std::istringstream iss(instr);
    float r,g,b,a;
    if( !(iss >> r >> g >> b >> a)){
      std::cout<<"Error reading transfer function data!"<<std::endl;
      break;
    }
    this->colortf[i] = cyColor(r,g,b);
    this->alphatf[i] = a;
    i++;
  }
  infile.close();

  return true;
  
  
}
/*Using this to parse data and find range; min; max */
void VolumeData::CreateHistogram(uchar* volume_data, int size)
{
  unsigned char min = 255;
  unsigned char max = 0;
  
  for(int i=0; i<size; i++){
    if(  volume_data[i] < min ) min = volume_data[i];
    if(  volume_data[i] > max ) max = volume_data[i];
  }
  std::cout<< "Min: "<<(uint)min<<"\nMax: "<<(uint)max<<std::endl;
  this->minData = min;
  this->maxData = max;
  return;
}

void VolumeData::CreateTransferFunction()
{
  colortf = new cyColor[this->n_bins];
  alphatf = new float[this->n_bins];

  /* assigning a single color and alpha value to each bin for now; need a better transfer function*/
  for(int i = 0; i < n_bins; i++){
    if ( i <= (int)n_bins/2 ){
	colortf[i] = cyColor(125, 35, 0);
	alphatf[i] = 0.01;
	if(i >= 15 && i <= 17) {
	  colortf[i] = cyColor(12,150,12);
	  alphatf[i] = 0.4;
	}
      }
      else{
	colortf[i] = cyColor(0, 160, 189);
	alphatf[i] = 0.0;
      }
  }
}

void VolumeData::GetTransferFunction(cyColor** color_tf, float** alpha_tf, int &tf_size, unsigned char &min, unsigned char &max)
{
  *color_tf = this->colortf;
  *alpha_tf = this->alphatf;
  tf_size = this->n_bins;
  min = this->minData;
  max = this->maxData;
}
