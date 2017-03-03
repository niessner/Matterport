
#pragma once

#include "mLibInclude.h"

class NormalExtractor {
public:

	//normals from bilateral filt depth
	static void computeDepthNormals(const SensorData* sd, const float depthSigmaD, const float depthSigmaR,
		unsigned int outWidth, unsigned int outHeight, std::vector<PointImage>& normals);

	static void computeMeshNormals(const SensorData* sd, const MeshDataf& mesh, 
		unsigned int outWidth, unsigned int outHeight, std::vector<PointImage>& normals) {
		throw MLIB_EXCEPTION("unimplemented");
	}

	static void saveNormalImage(const std::string& filename, const PointImage& normal)
	{
		//convert to 3x 16bit
		BaseImage<vec3us> normalImage(normal.getDimensions());
		for (auto& p : normal) {
			vec3f n = (p.value.x == -std::numeric_limits<float>::infinity()) ? vec3f(0.0f, 0.0f, 0.0f) : p.value;
			n = 0.5f * (n + 1) * 65535.0f; //[-1,1] -> [0,65535]
			normalImage(p.x, p.y) = vec3us(n);
		}

		FreeImageWrapper::saveImage(filename, normalImage);
	}

private:

};
