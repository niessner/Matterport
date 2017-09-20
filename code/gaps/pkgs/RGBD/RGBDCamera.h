// Include file for RGBDCamera class



////////////////////////////////////////////////////////////////////////
// Class definition
////////////////////////////////////////////////////////////////////////

class RGBDCamera {
public:
  // Constructors/destructors
  RGBDCamera(void);
  RGBDCamera(const RGBDCamera& camera);
  ~RGBDCamera(void);

  // Parameter manipulation functions
  void SetUndistortedParameters(void);

  // Input/output functions
  int ReadScanNetFile(const char *filename);
  int ReadMatterportFile(const char *filename);
  int WriteScanNetFile(const char *filename) const;
  int WriteMatterportFile(const char *filename) const;

public:
  // Input data
  int colorWidth;    
  int colorHeight; 
  int depthWidth;  
  int depthHeight; 
  double fx_color;
  double fy_color;
  double mx_color; 
  double my_color; 
  double k1_color;
  double k2_color;
  double k3_color;
  double k4_color;
  double p1_color;
  double p2_color;
  double fx_depth; 
  double fy_depth; 
  double mx_depth; 
  double my_depth;
  double k1_depth;
  double k2_depth;
  double k3_depth;
  double k4_depth;
  double p1_depth;
  double p2_depth;
  R4Matrix depthToColorExtrinsics;
  char *deviceId;          
  char *deviceName;        
  char *sceneLabel;        
  char *sceneType;         
  int numDepthFrames;      
  int numColorFrames;      
  int numIMUmeasurements;  
};



