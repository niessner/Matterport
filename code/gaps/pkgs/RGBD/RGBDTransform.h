////////////////////////////////////////////////////////////////////////
// Include file for RGBD coordinate transformation functions
////////////////////////////////////////////////////////////////////////

extern int RGBDTransformTextureToSurface(const R2Point& texture_position, R2Point& surface_position, const RGBDSurface *surface);
extern int RGBDTransformTextureToWorld(const R2Point& texture_position, R3Point& world_position, const RGBDSurface *surface);
extern int RGBDTransformTextureToCamera(const R2Point& texture_position, R3Point& camera_position, const RGBDSurface *surface, const RGBDImage *image);
extern int RGBDTransformTextureToImage(const R2Point& texture_position, R2Point& image_position, const RGBDSurface *surface, const RGBDImage *image);

extern int RGBDTransformSurfaceToTexture(const R2Point& surface_position, R2Point& texture_position, const RGBDSurface *surface);
extern int RGBDTransformSurfaceToWorld(const R2Point& surface_position, R3Point& world_position, const RGBDSurface *surface);
extern int RGBDTransformSurfaceToCamera(const R2Point& surface_position, R3Point& camera_position, const RGBDSurface *surface, const RGBDImage *image);
extern int RGBDTransformSurfaceToImage(const R2Point& surface_position, R2Point& image_position, const RGBDSurface *surface, const RGBDImage *image);

extern int RGBDTransformWorldToTexture(const R3Point& world_position, R2Point& texture_position, const RGBDSurface *surface);
extern int RGBDTransformWorldToSurface(const R3Point& world_position, R2Point& surface_position, const RGBDSurface *surface);
extern int RGBDTransformWorldToCamera(const R3Point& world_position, R3Point& camera_position, const RGBDImage *image);
extern int RGBDTransformWorldToImage(const R3Point& world_position, R2Point& image_position, const RGBDImage *image);

extern int RGBDTransformCameraToTexture(const R3Point& camera_position, R2Point& texture_position, const RGBDImage *image, const RGBDSurface *surface);
extern int RGBDTransformCameraToSurface(const R3Point& camera_position, R2Point& surface_position, const RGBDImage *image, const RGBDSurface *surface);
extern int RGBDTransformCameraToWorld(const R3Point& camera_position, R3Point& wold_position, const RGBDImage *image);
extern int RGBDTransformCameraToImage(const R3Point& camera_position, R2Point& image_position, const RGBDImage *image);

extern int RGBDTransformImageToTexture(const R2Point& image_position, R2Point& texture_position, const RGBDImage *image, const RGBDSurface *surface);
extern int RGBDTransformImageToSurface(const R2Point& image_position, R2Point& surface_position, const RGBDImage *image, const RGBDSurface *surface);
extern int RGBDTransformImageToWorld(const R2Point& image_position, R3Point& world_position, const RGBDImage *image);
extern int RGBDTransformImageToCamera(const R2Point& image_position, R3Point& camera_position, const RGBDImage *image);

