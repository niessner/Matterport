
#pragma once

#include "mLibInclude.h"

#define FIND_FAR_CAMERAS

struct ReprojError {
	size_t numCorrs;
	float depthL1;
	float depthL2;
	float intensityL1;
	float intensityGradL1;

	ReprojError() {
		numCorrs = 0;
		depthL1 = 0;
		depthL2 = 0;
		intensityL1 = 0;
		intensityGradL1 = 0;
	}

	ReprojError operator+(const ReprojError& other) const {
		ReprojError ret;
		ret.numCorrs = numCorrs + other.numCorrs;
		ret.depthL1 = depthL1 + other.depthL1;
		ret.depthL2 = depthL2 + other.depthL2;
		ret.intensityL1 = intensityL1 + other.intensityL1;
		ret.intensityGradL1 = intensityGradL1 + other.intensityGradL1;
		return ret;
	}
	inline void operator+=(const ReprojError& other) {
		numCorrs += other.numCorrs;
		depthL1 += other.depthL1;
		depthL2 += other.depthL2;
		intensityL1 += other.intensityL1;
		intensityGradL1 += other.intensityGradL1;
	}

	void normalize() {
		if (numCorrs > 0) {
			depthL1 /= (float)numCorrs;
			depthL2 /= (float)numCorrs;
			intensityL1 /= (float)numCorrs;
			intensityGradL1 /= (float)numCorrs;
		}
	}
};


class ScannedScene {
public:
	ScannedScene(const std::string& path, const std::string& name) {
		load(path, name);
	}
	~ScannedScene() {
		for (auto* sd : m_sds) {
			SAFE_DELETE(sd);
		}
	}

	void load(const std::string& path, const std::string& name) {

		m_name = name;
		m_path = path;

		Directory dir(path);
		auto& files = dir.getFilesWithSuffix(".sens");
		std::sort(files.begin(), files.end());

		std::cout << name << std::endl;
		for (size_t i = 0; i < files.size(); i++) {
			auto& f = files[i];
			m_sds.push_back(new SensorData);
			SensorData* sd = m_sds.back();
			sd->loadFromFile(path + "/" + f);
			std::cout << *sd << std::endl;
		}
	}

	ReprojError computeReprojection(size_t maxNumFramePairSamples, size_t maxNumSampleTries, float maxDepth, std::vector<float>& sampledCameraDists);


private:
	std::vector<SensorData*> m_sds;
	std::string m_name;
	std::string m_path;

};

class Reprojection {
public:
	static ReprojError computeReprojection(const DepthImage32& depth0, const ColorImageR32& intensity0,
		const DepthImage32& depth1, const ColorImageR32& intensity1, const mat4f& depthIntrinsics0, const mat4f& depthIntrinsics1,
		const mat4f& depthIntrinsicsInv0, const mat4f& depthIntrinsicsInv1, const mat4f& transform0to1, float maxDepth);

	static vec3f computeNormal(const DepthImage32& depth, const mat4f& intrinsicsInv, unsigned int x, unsigned int y)
	{
		if (x > 0 && x + 1 < depth.getWidth() && y > 0 && y + 1 < depth.getHeight()) {
			const vec3f& CC = computeCameraSpacePosition(depth, intrinsicsInv, x + 0, y + 0);
			const vec3f& PC = computeCameraSpacePosition(depth, intrinsicsInv, x + 0, y + 1);
			const vec3f& CP = computeCameraSpacePosition(depth, intrinsicsInv, x + 1, y + 0);
			const vec3f& MC = computeCameraSpacePosition(depth, intrinsicsInv, x + 0, y - 1);
			const vec3f& CM = computeCameraSpacePosition(depth, intrinsicsInv, x - 1, y + 0);

			if (CC.x != -std::numeric_limits<float>::infinity() && PC.x != -std::numeric_limits<float>::infinity() &&
				CP.x != -std::numeric_limits<float>::infinity() && MC.x != -std::numeric_limits<float>::infinity() &&
				CM.x != -std::numeric_limits<float>::infinity())
			{
				const vec3f n = (PC - MC) ^ (CP - CM);
				const float l = n.length();
				if (l > 0.0f) return n / -l;
				else return vec3f(-std::numeric_limits<float>::infinity());
			}
		}
		return vec3f(-std::numeric_limits<float>::infinity());
	}

	static vec3f computeCameraSpacePosition(const DepthImage32& depth, const mat4f& intrinsicsInv, unsigned int x, unsigned int y)
	{
		const float d = depth(x, y);
		if (d == -std::numeric_limits<float>::infinity() || d == 0) return vec3f(-std::numeric_limits<float>::infinity());
		return (intrinsicsInv*vec4f((float)x*d, (float)y*d, d, d)).getVec3();
	}

	static vec2f cameraToDepth(const mat4f& depthIntrinsics, const vec3f& pos)
	{
		vec3f p = depthIntrinsics * pos;
		return vec2f(p.x / p.z, p.y / p.z);
	}

	//static vec3f cameraToKinectProj(const vec3f& pos, const mat4f& intrinsics, unsigned int width, unsigned int height, float maxDepth)	{
	//	vec2f proj = cameraToDepth(intrinsics, pos);
	//	vec3f pImage = vec3f(proj.x, proj.y, pos.z);
	//	pImage.x = (2.0f*pImage.x - (width - 1.0f)) / (height - 1.0f);
	//	pImage.y = ((height - 1.0f) - 2.0f*pImage.y) / (height - 1.0f);
	//	pImage.z = (pImage.z - 0.4f) / (maxDepth - 0.4f);
	//	return pImage;
	//}

	static bool hasCameraFrustumIntersection(const mat4f& t0, const mat4f& intrinsics0, const mat4f& t1, const mat4f& intrinsics1,
		unsigned int width, unsigned int height, float maxDepth)
	{
		const float cameraDist = vec3f::dist(t0.getTranslation(), t1.getTranslation());
		if (cameraDist > 6.0f) //maxDepth) 
			return false;
		const vec3f d0 = t0.getRotation() * vec3f(1.0f, 1.0f, 1.0f).getNormalized();
		const vec3f d1 = t1.getRotation() * vec3f(1.0f, 1.0f, 1.0f).getNormalized();
		const float cameraAngles = std::acos(math::clamp(d0 | d1, -1.0f, 1.0f));
		if (cameraAngles > 2.0f)
			return false;

		const float minDepth = 0.4f;

		const vec2f imgMin(0.0f); const vec2f imgMax((float)width, (float)height);
		std::vector<vec3f> corners(4);
		corners[0] = vec3f(imgMin.x, imgMin.y, 1.0f);		corners[1] = vec3f(imgMin.x, imgMax.y, 1.0f);
		corners[2] = vec3f(imgMax.x, imgMax.y, 1.0f);		corners[3] = vec3f(imgMax.x, imgMin.y, 1.0f);

		mat4f t0_inv = t0.getInverse();
		mat4f t1_inv = t1.getInverse();
		const mat4f intrinsicsInv0 = intrinsics0.getInverse();
		const mat4f intrinsicsInv1 = intrinsics1.getInverse();

		BoundingBox2f i_in_j, j_in_i;
		for (unsigned int c = 0; c < corners.size(); c++) {
			vec2f ij = cameraToDepth(intrinsics1, t1_inv * t0 * intrinsicsInv0 * (minDepth * corners[c]));
			vec2f ji = cameraToDepth(intrinsics0, t0_inv * t1 * intrinsicsInv1 * (minDepth * corners[c]));
			i_in_j.include(ij);
			j_in_i.include(ji);

			ij = cameraToDepth(intrinsics1, t1_inv * t0 * intrinsicsInv0 * (4.0f * corners[c]));//(maxDepth * corners[c]));
			ji = cameraToDepth(intrinsics0, t0_inv * t1 * intrinsicsInv1 * (4.0f * corners[c]));//(maxDepth * corners[c]));
			i_in_j.include(ij);
			j_in_i.include(ji);
		}
		BoundingBox2f bbImg(imgMin, imgMax);
		if (!(bbImg.intersects(i_in_j) && bbImg.intersects(j_in_i)))
			return false;

		return true;
	}

	static bool checkValidReprojection(const DepthImage32& depth0, const DepthImage32& depth1,
		const mat4f& depthIntrinsics0, const mat4f& depthIntrinsics1,
		const mat4f& depthIntrinsicsInv0, const mat4f& depthIntrinsicsInv1, const mat4f& transform0to1, float maxDepth);
private:
};