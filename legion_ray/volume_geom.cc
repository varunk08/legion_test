#include "volume_geom.h"
//-----------------------------------------------------------------------------------------------------

  VolumeGeometry::VolumeGeometry(int x_width, int y_width, int z_width,
				 const std::vector<PhysicalRegion> &regions, Context ctx, HighLevelRuntime *runtime)
  {
    pmin = Point3(-1, -1, -1);
    pmax = Point3(1, 1, 1);
    this->data_xdim = x_width;
    this->data_ydim = y_width;
    this->data_zdim = z_width;
    this->vol_regions = regions;
    this->ctx = ctx;
    this->runtime = runtime;
  }
//-----------------------------------------------------------------------------------------------------
  Box VolumeGeometry::GetBoundBox() const
  {
    return Box(-1,-1,-1,1,1,1);
  }
//-----------------------------------------------------------------------------------------------------
void VolumeGeometry::init_logical_regions(LogicalRegion vol_data_lr, LogicalRegion tf_lr,
					  GucharRA &vol_acc, GPoint3fRA &acc_grads,
					  GRGBColorRA &acc_col, GfloatRA &acc_alpha)
  {

    this->acc_vol_data = vol_acc;
    this->acc_gradients = acc_grads;
    this->acc_color_tf = acc_col;
    this->acc_alpha_tf = acc_alpha;

    this->corner_pos_lr = vol_data_lr;
    this->tf_lr = tf_lr;

    
  }
//-----------------------------------------------------------------------------------------------------
void VolumeGeometry::init_tf_bounds(int min, int max)
{
  this->minData = min;
  this->maxData = max;
}
//-----------------------------------------------------------------------------------------------------
void VolumeGeometry::set_lights(LightList &lList)
{
 PointLight *light = new PointLight();
  AmbientLight *amb_light = new AmbientLight();
 
 
  amb_light->SetIntensity(cyColor(0.2,0.2,0.2));
  light->SetIntensity(cyColor(0.7,0.7,0.7));
  light->SetPosition(Point3(0,10,0));
  this->lights.push_back(amb_light);
  this->lights.push_back(light);
 
  /*  this->light = new PointLight;
  this->light->SetIntensity(cyColor(0.7,0.7,0.7));
  this->light->SetPosition(Point3(0,10,0));*/
 
}
//-----------------------------------------------------------------------------------------------------
  void VolumeGeometry::SetTransform(const Matrix3 &nodeToWorld, const Matrix3 &n2w_itm, const Point3 &pos ){
    this->tm.Set(nodeToWorld.data);
    this->itm.Set(n2w_itm.data);
    this->worldPos = Point3(pos);
  }
//-----------------------------------------------------------------------------------------------------
  bool VolumeGeometry::IntersectRay(const Ray &ray, HitInfo &hInfo, int hitSide) const
  {
    bool boxHit = false;
    //the ray is already in object space.
    Point3 corners[8];
    corners[0] = Point3(-1.0f, -1.0f, -1.0f);
    corners[1] = Point3(1.0f, -1.0f, -1.0f);
    corners[2] = Point3(1.0f, -1.0f, 1.0f);
    corners[3] = Point3(-1, -1, 1);
    corners[4] = Point3(-1, 1, -1);
    corners[5] = Point3(1, 1, -1);
    corners[6] = Point3(1, 1, 1);
    corners[7] = Point3(-1, 1, 1);
    int indices[36] = {
      0,1,2,
      0,2,3,//bottom
      3,0,4,
      3,4,7,//left
      4,5,6,
      4,6,7,//top
      5,2,6,
      5,1,2,//right
      0,5,4,
      0,1,5,//front
      2,3,7,
      2,7,6//back
    };
    HitInfo tempHInfo;
    tempHInfo.Init();
    tempHInfo.node = hInfo.node;
    tempHInfo.z = hInfo.z;
    std::vector<HitInfo> hitZ(2);
    hitZ.clear();

    for(int i = 0; i < 36; i+=3){
      if( Box_IntersectTriangle ( ray, tempHInfo, corners[indices[i]], corners[indices[i+1]], corners[indices[i+2]] )){
	hitZ.push_back(tempHInfo);
	boxHit = true;
	tempHInfo.Init();
	tempHInfo.node = hInfo.node;
	tempHInfo.z = hInfo.z;
      }
    }

    if(boxHit) {
      int st = 0, end = 0;
      if( hitZ[0].z < hInfo.z && hitZ[1].z < hInfo.z ) {
	if( hitZ[0].z < hitZ[1].z ) { 
	  st = 0; end = 1;
	}
	else {
	  st = 1; end = 0;
	}
      }
      bool noData = false;
      hInfo.shade = RayMarch(ray, hInfo, hitZ[st], hitZ[end], noData);
      if(!noData) {
	hInfo.front =  true;
      }
      else boxHit = false; //no data (or) iso surface rendering
    }

    return boxHit;
  }

  //-----------------------------------------------------------------------------------------------------

  cyColor VolumeGeometry::RayMarch(const Ray& ray, HitInfo &hInfo, HitInfo& start, HitInfo& end, bool &noData) const
  {

    cyColor shade = cyColor(0,0,0); 
    float length = (float)end.z - start.z;
    float dist = sqrt(std::pow((float)end.p.x - start.p.x, 2.0f) + pow((float)end.p.y - start.p.y, 2.0f) + pow((float)end.p.z - start.p.z,  2.0f));

    /** Binning the sample point */
    Ray r;
    r.p = start.p;
    r.dir = end.p - start.p;
    r.Normalize();
    float dt = (float)dist / (data_xdim * 2.0f ); 
    float t = 0;
    Point3 sample;
    cyColor color_acc = cyColor(0,0,0);
    float alpha_acc = 0.0f;
    int num_iter = 0;
    noData = true;

    //on the fly calculation of corner positions
    float startX = -1.0f;
    float startY = -1.0f;
    float startZ = -1.0f;

    float dx = (2.0f / (data_xdim - 1));
    float dy = (2.0f / (data_ydim - 1));
    float dz = (2.0f / (data_zdim - 1));

    while ( t < dist ) {
      sample = start.p + t * r.dir;
      int cx=-1,cy=-1,cz=-1;

      float low_bnd;
      float up_bnd;

      for(int x=0; x < data_xdim -1; x++){
	low_bnd = startX + x * dx;
	up_bnd = startX + (x+1) * dx;

	if((float)low_bnd < sample.x && (float)up_bnd >= sample.x){
	  cx = x;
	  break;
	}

      }

      for(int y=0; y < data_ydim-1; y++){
	low_bnd = startY + y * dy;
	up_bnd = startY + (y+1) * dy;
	if(low_bnd < sample.y && up_bnd >= sample.y){
	  cy = y;
	  break;
	}
      }
     
      for(int z=0; z < data_zdim-1; z++){
	low_bnd = startZ + z * dz;
	up_bnd = startZ + (z+1) * dz;

	if(low_bnd < sample.z && up_bnd >= sample.z){
	  cz = z;	  
	  break;
	}

      }

      if(cx >=0 && cy >=0 && cz >=0 ){
	num_iter++;
	float data_tot = TrilinearInterpolate(sample,cx,cy,cz);
	bool compositing =true;
	if( compositing ){
	  Point3 P = sample;
	  Point3 N = EstimateGradient(cx, cy, cz);
	  cyColor cin = DoPhongShading(ray, P, N, GetColor(data_tot));
	  float ain = GetAlpha(data_tot);
	  //color accumulation; cacc += (1-opacc)*color*opacity;
	  //Do phong shading
	  color_acc += (1 - alpha_acc) * cin * ain;
	  //alpha accumulation; opacc += (1- opacc)*opacity;
	  alpha_acc += (1 - alpha_acc) * ain;
	
	  //compositing
	  if( alpha_acc > 0.9 ){
	    hInfo.p = P;
	    hInfo.N = N;
	    noData = false;
	    hInfo.z += t;
	    hInfo.volume = true;
	    shade = color_acc;
	    break;
	  }
	}//compositing
	else {
	  //iso surface rendering
	  float THRESHOLD = 80;
	  if( data_tot > THRESHOLD ) {
	    noData = false;
	    hInfo.p = sample;
	    hInfo.N = EstimateGradient(cx,cy,cz);
	    hInfo.z += t;
	    hInfo.renderIsoSurface = true;
	    hInfo.volume = false;
	    break;
	  }
	}//iso surface rendering
      }//valid indices cx, cy, cz

      t += dt;
    }//while loop
    
    return shade;
    
  }
  //-----------------------------------------------------------------------------------------------------
  /*Shading needs: ray, lights, colors(amb, diff, spec)*/
  /*Lights are in world coords!! */
cyColor VolumeGeometry::DoPhongShading(const Ray &ray,const Point3 &pt, const Point3 &norm, cyColor preCol) const
  {
    /*
      Increase ambient
      Fix specular
    */
    cyColor shade;
    float amp = 3.0f;
    cyColor diffuse =  amp * preCol;
    cyColor ambInt =  amp * preCol;
    cyColor ambComponent = cyColor(0,0,0);
    cyColor specular = cyColor(0.7f,0.7f,0.7f);
    cyColor allOther = cyColor(0,0,0);
    float glossiness = 20.0f;    
    Point3 P = TransformFrom(pt);
    Point3 N = VectorTransformFrom(norm);
    N.Normalize();
  
  for ( unsigned int i=0; i<this->lights.size(); i++ ) {
      if(lights[i]->IsAmbient()){
	cyColor intensity = this->lights[i]->Illuminate(P);
	ambComponent += (ambInt * intensity);
	continue;
      }
      else{
    
	Point3 L = -lights[i]->Direction(P);
	L.Normalize();
            
	Point3 V = TransformFrom(ray.p) - P;
	V.Normalize();
            
	Point3 LplusV = L + V;
	Point3 H = (L+V)/LplusV.Length();
	H.Normalize();

	float costheta = L.Dot(N)/(L.Length() * N.Length());
	float alpha = glossiness;

	float S = N.Dot(H);
	S = (S > 0)?S:0.0f; //max
	S = pow( (float)S , alpha );

	cyColor intensity = lights[i]->Illuminate(P);
	allOther += intensity *  ((costheta>0?costheta:0)  * diffuse) + (S * specular);
      }
      // finally add inta*cola + intall*costheta*(cold + s* colS)
      shade = ambComponent  + allOther;
    }
    return shade;
  }
  //-----------------------------------------------------------------------------------------------------
  cyColor VolumeGeometry::GetColor(float density) const
  {
    RGBColor col;
    float f_index = (density - (unsigned int)minData)/((unsigned int)maxData-(unsigned int)minData) * 255;
    if(f_index < 255){
      int first_ind = (int) f_index;
      int next_ind = (int) 1 + first_ind;
      float alpha = f_index - first_ind;

    RGBColor v0 = acc_color_tf.read(DomainPoint::from_point<1>(Point<1>(first_ind)));
    RGBColor v1 = acc_color_tf.read(DomainPoint::from_point<1>(Point<1>(next_ind)));

    col.r = (1 - alpha) * v0.r + alpha * v1.r;
    col.g = (1 - alpha) * v0.g + alpha * v1.g;
    col.b = (1 - alpha) * v0.b + alpha * v1.b;

    }
    else{
      col = acc_color_tf.read(DomainPoint::from_point<1>(Point<1>( (int) f_index) ));
    }

   
    return cyColor(col.r, col.g, col.b);
  }
  //-----------------------------------------------------------------------------------------------------
  float VolumeGeometry::GetAlpha(float density) const
  {
    int index;
    float lerp_alpha;
    float f_index = (density - (unsigned int)minData)/((unsigned int)maxData-(unsigned int)minData) * 255;

    if(f_index < 255){
      int first_ind = (int) f_index;
      int next_ind = (int) 1 + first_ind;
      float alpha = f_index - first_ind;

    float v0 = acc_alpha_tf.read(DomainPoint::from_point<1>(Point<1>(first_ind)));
    float v1 = acc_alpha_tf.read(DomainPoint::from_point<1>(Point<1>(next_ind)));

    lerp_alpha = (1 - alpha) * v0 + alpha * v1;


    }
    else{
      lerp_alpha = acc_alpha_tf.read(DomainPoint::from_point<1>(Point<1>( (int) f_index) ));
    }

    return lerp_alpha;
  }
  //-----------------------------------------------------------------------------------------------------
  inline int VolumeGeometry::getIndex(int x, int y, int z) const
  {
    return ( z * this->data_ydim + y ) * this->data_xdim + x;
  }
  //-----------------------------------------------------------------------------------------------------

  float VolumeGeometry::TrilinearInterpolate(Point3 samplePt, int x, int y, int z) const
  {
    float data = 0;
    float xd, yd, zd;
    float c01, c23, c45, c67;
    float c0,c1;
    float val;

    float startX = -1.0f;
    float startY = -1.0f;
    float startZ = -1.0f;

    float dx = (2.0f / (data_xdim - 1));
    float dy = (2.0f / (data_ydim - 1));
    float dz = (2.0f / (data_zdim - 1));


    Point3f t_p1;
    t_p1.x = startX + x * dx;
    t_p1.y = startY + y * dy;
    t_p1.z = startZ + z * dz; 
    Point3f t_p2;
    t_p2.x = startX + (x+1) * dx;
    t_p2.y = startY + (y+1) * dy;
    t_p2.z = startZ + (z+1) * dz;
    float* v = new float[8];

    v[0] = (unsigned int)acc_vol_data.read(DomainPoint::from_point<1>(Point<1>(getIndex(x,y,z))));
    v[1] = (unsigned int)acc_vol_data.read(DomainPoint::from_point<1>(Point<1>(getIndex(x+1,y,z))));
    v[2] = (unsigned int)acc_vol_data.read(DomainPoint::from_point<1>(Point<1>(getIndex(x,y+1,z))));
    v[3] = (unsigned int)acc_vol_data.read(DomainPoint::from_point<1>(Point<1>(getIndex(x+1,y+1,z))));
    v[4] = (unsigned int)acc_vol_data.read(DomainPoint::from_point<1>(Point<1>(getIndex(x+1,y,z+1))));
    v[5] = (unsigned int)acc_vol_data.read(DomainPoint::from_point<1>(Point<1>(getIndex(x,y,z+1))));
    v[6] = (unsigned int)acc_vol_data.read(DomainPoint::from_point<1>(Point<1>(getIndex(x,y+1,z+1))));
    v[7] = (unsigned int)acc_vol_data.read(DomainPoint::from_point<1>(Point<1>(getIndex(x+1,y+1,z+1))));

    Point3 p1, p2;
    p1.x = t_p1.x;     p1.y = t_p1.y;     p1.z = t_p1.z; 
    p2.x = t_p2.x;     p2.y = t_p2.y;     p2.z = t_p2.z; 
    xd = (samplePt.x - p1.x)/(p2.x - p1.x);
    yd = (samplePt.y - p1.y)/(p2.y - p1.y);
    zd = (samplePt.z - p1.z)/(p2.z - p1.z);
    c01 = v[0] * (1 - xd) + v[1] * xd;
    c23 = v[2] * (1 - xd) + v[3] * xd;
    c45 = v[4] * (1 - xd) + v[5] * xd;
    c67 = v[6] * (1 - xd) + v[7] * xd;
    c0 = c01 * (1 - yd) + c45 * yd;
    c1 = c23 * (1 - yd) + c67 * yd;
    val = c0 * (1 - zd) + c1 * zd;

    delete[] v;
    return val;
  }
  //-----------------------------------------------------------------------------------------------------
  Point3 VolumeGeometry::EstimateGradient(int x, int y, int z) const
  {

    /* maybe i can calculate this avg while calculating gradients; maybe needto trilinear interpolate gradient
       as well */
    Point3 grad(0,0,0);
    Point3f acc_g;
    Point3f temp =  (Point3f)acc_gradients.read(DomainPoint::from_point<1>(Point<1>(getIndex(x,y,z))));
    acc_g.x = temp.x;  acc_g.y = temp.y;  acc_g.z = temp.z; 
    temp=  acc_gradients.read(DomainPoint::from_point<1>(Point<1>(getIndex(x+1,y,z))));
    acc_g.x += temp.x;  acc_g.y += temp.y;  acc_g.z += temp.z; 
    temp=  acc_gradients.read(DomainPoint::from_point<1>(Point<1>(getIndex(x,y+1,z))));
    acc_g.x += temp.x;  acc_g.y += temp.y;  acc_g.z += temp.z; 
    temp=  acc_gradients.read(DomainPoint::from_point<1>(Point<1>(getIndex(x+1,y+1,z))));
    acc_g.x += temp.x;  acc_g.y += temp.y;  acc_g.z += temp.z; 
    temp=  acc_gradients.read(DomainPoint::from_point<1>(Point<1>(getIndex(x+1,y,z+1))));
    acc_g.x += temp.x;  acc_g.y += temp.y;  acc_g.z += temp.z; 
    temp=  acc_gradients.read(DomainPoint::from_point<1>(Point<1>(getIndex(x,y,z+1))));
    acc_g.x += temp.x;  acc_g.y += temp.y;  acc_g.z += temp.z; 
    temp=  acc_gradients.read(DomainPoint::from_point<1>(Point<1>(getIndex(x,y+1,z+1))));
    acc_g.x += temp.x;  acc_g.y += temp.y;  acc_g.z += temp.z; 
    temp=  acc_gradients.read(DomainPoint::from_point<1>(Point<1>(getIndex(x+1,y+1,z+1))));
    acc_g.x += temp.x;  acc_g.y += temp.y;  acc_g.z += temp.z; 
  
    grad.x = acc_g.x;    grad.y = acc_g.y;    grad.z = acc_g.z;
    grad.x /= -8.0; grad.y /= -8.0; grad.z /= -8.0;
    if(grad.x != 0 && grad.y!= 0 && grad.z!=0)
      grad.Normalize();

    return grad;
    
  }

  //-----------------------------------------------------------------------------------------------------
  bool VolumeGeometry::Box_IntersectTriangle(const Ray &ray, HitInfo &hInfo, Point3 p1, Point3 p2, Point3 p3, bool raymarch) const
  {
    //Based on Shirley's book
    Point3 A, B, C;
    A = p1; B = p2; C = p3;
    float a = A.x - B.x;
    float b = A.y - B.y;
    float c = A.z - B.z;

    float d = A.x - C.x;
    float e = A.y - C.y;
    float f = A.z - C.z;

    float g = ray.dir.x;
    float h = ray.dir.y;
    float i = ray.dir.z;

    float j = A.x - ray.p.x;
    float k = A.y - ray.p.y;
    float l = A.z - ray.p.z;

    float eimhf = e * i - h * f;
    float gfmdi = g * f - d * i;
    float dhmeg = d * h - e * g;
    float akmjb = a * k - j * b;
    float jcmal = j * c - a * l;
    float blmkc = b * l - k * c;

    float M = a * (eimhf) + b * (gfmdi) + c * (dhmeg);
    float t = -(f * akmjb + e * jcmal + d * blmkc)/M;

    if( t > BIGFLOAT || t > hInfo.z)
      return false;

    float gamma = (i*(akmjb) + h * ( jcmal ) + g * ( blmkc ))/M;

    if( gamma  < 0 || gamma > 1)
      return false;

    float beta = (j * (eimhf) + k * gfmdi + l * dhmeg )/M;

    if(beta < 0 || beta > 1 - gamma)
      return false;

    hInfo.z = t;
    hInfo.p = ray.p + t * ray.dir;
    Point3 N = (C - A).Cross(B - A);
    N.Normalize();
    hInfo.N = N;
    hInfo.volume = true;

    
    return true;
  }
//-----------------------------------------------------------------------------------------------------
  Point3 VolumeGeometry::TransposeMult( const Matrix3 &m, const Point3 &dir ) const
  {
    Point3 d;
    d.x = m.GetColumn(0) % dir;
    d.y = m.GetColumn(1) % dir;
    d.z = m.GetColumn(2) % dir;
    return d;
  }

//-----------------------------------------------------------------------------------------------------

// Transform from the local coordinate system
Point3 VolumeGeometry::TransformFrom( const Point3 &p ) const { return (tm) * (p) + worldPos; }

//-----------------------------------------------------------------------------------------------------

Point3 VolumeGeometry::VectorTransformFrom( const Point3 &dir ) const { return TransposeMult(itm,dir); }

//-----------------------------------------------------------------------------------------------------
