
#include "stdafx.h"

#include "normals.h"
#include "imageHelper.h"


void NormalExtractor::computeDepthNormals(const SensorData* sd, float depthSigmaD, float depthSigmaR, std::vector<PointImage>& normals)
{
	normals.resize(sd->m_frames.size());

	mat4f intrinsic = sd->m_calibrationDepth.m_intrinsic;
	if (m_width != sd->m_depthWidth || m_height != sd->m_depthHeight) {
		// adapt depth intrinsics
		intrinsic._m00 *= (float)m_width / (float)sd->m_depthWidth;				//focal length
		intrinsic._m11 *= (float)m_height / (float)sd->m_depthHeight;			//focal length
		intrinsic._m02 *= (float)(m_width - 1) / (float)(sd->m_depthWidth - 1);		//principal point
		intrinsic._m12 *= (float)(m_height - 1) / (float)(sd->m_depthHeight - 1);	//principal point
	}
	const mat4f depthIntrinsicInv = intrinsic.getInverse();

#pragma omp parallel for
	for (int imageIdx = 0; imageIdx < (int)sd->m_frames.size(); imageIdx++) {
		DepthImage32 d = sd->computeDepthImage(imageIdx);
		//DepthImage32 dfilt = d; 
		//ImageHelper::bilateralFilter(dfilt, depthSigmaD, depthSigmaR);
		//PointImage normal = SensorData::computeNormals(depthIntrinsicInv, dfilt);
		//normal.resize(m_width, m_height);
		//normals[imageIdx] = normal;

		d.resize(m_width, m_height);
		DepthImage32 dfilt = d;
		ImageHelper::bilateralFilter(dfilt, depthSigmaD, depthSigmaR);
		normals[imageIdx] = SensorData::computeNormals(depthIntrinsicInv, dfilt);

		////debugging
		//std::string file = "t.png";
		//saveNormalImage(file, normals[imageIdx]);
		//PointImage check;
		//NormalExtractor::loadNormalImage(file, check);
		//for (const auto& p : check) {
		//	const auto& n = normals[imageIdx](p.x, p.y);
		//	if (p.value.x == -std::numeric_limits<float>::infinity() && n.x != -std::numeric_limits<float>::infinity())
		//		int  a = 5;
		//	else if (p.value.x != -std::numeric_limits<float>::infinity() && n.x == -std::numeric_limits<float>::infinity())
		//		int  a = 5;
		//	else {
		//		float d = vec3f::dist(p.value, n);
		//		if (d > 0.0001f)
		//			int a = 5;
		//	}
		//}
		//PointImage campos = sd->computeCameraSpacePositions((unsigned int)imageIdx);
		//campos.resize(m_width, m_height);
		//{
		//	PointCloudf pc;
		//	for (const auto& p : campos) {
		//		const vec3f& n = check(p.x, p.y);
		//		if (p.value.x != -std::numeric_limits<float>::infinity() && n.x != -std::numeric_limits<float>::infinity()) {
		//			pc.m_points.push_back(p.value);
		//			pc.m_normals.push_back(n);
		//		}
		//	}
		//	PointCloudIOf::saveToFile("t.ply", pc);
		//}
		//int a = 5;
		////debugging
	}
}

void NormalExtractor::computeMeshNormals(const SensorData* sd, const TriMeshf& mesh, float zNear, float zFar, std::vector<PointImage>& normals)
{
	normals.resize(sd->m_frames.size());

	m_mesh.init(*m_graphics, mesh);

	mat4f intrinsic = sd->m_calibrationDepth.m_intrinsic;
	if (m_width != sd->m_depthWidth || m_height != sd->m_depthHeight) {
		// adapt depth intrinsics
		intrinsic._m00 *= (float)m_width / (float)sd->m_depthWidth;				//focal length
		intrinsic._m11 *= (float)m_height / (float)sd->m_depthHeight;			//focal length
		intrinsic._m02 *= (float)(m_width - 1) / (float)(sd->m_depthWidth - 1);		//principal point
		intrinsic._m12 *= (float)(m_height - 1) / (float)(sd->m_depthHeight - 1);	//principal point
	}

	for (size_t imageIdx = 0; imageIdx < sd->m_frames.size(); imageIdx++) {
		const mat4f& camToWorld = sd->m_frames[imageIdx].getCameraToWorld();
		const mat4f worldToCam = camToWorld.getInverse();
		const vec3f eye = camToWorld.getTranslation();
		const mat4f proj = projFromVision(intrinsic, camToWorld, m_width, m_height, zNear, zFar);

		ColorImageR32G32B32A32 image;
#pragma omp critical 
		{
			render(proj, worldToCam, eye, image, sd);
		}

		normals[imageIdx].allocate(image.getWidth(), image.getHeight());
		for (auto& p : normals[imageIdx]) {
			p.value = image(p.x, p.y).getVec3();
		}

		PointImage campos = sd->computeCameraSpacePositions((unsigned int)imageIdx);
		campos.resize(m_width, m_height);
		{
			PointCloudf pc;
			for (const auto& p : campos) {
				vec3f n = normals[imageIdx](p.x, p.y);
				if (p.value.x != -std::numeric_limits<float>::infinity() && (!(n.x == 0 && n.y == 0 && n.z == 0))) {
					n = 2.0f * n - 1.0f;
					if (std::abs(n.length() - 1) > 0.001f)
						int  a= 5;
					pc.m_points.push_back(p.value);
					pc.m_normals.push_back(n);
				}
			}
			PointCloudIOf::saveToFile("test.ply", pc);
		}

		FreeImageWrapper::saveImage("test.png", normals[imageIdx]);
		saveNormalImage("test2.png", normals[imageIdx]);
		int a = 5;
	}
}

ml::mat4f NormalExtractor::projFromVision(const mat4f& intrinsic, const mat4f& cameraToWorld, unsigned int width, unsigned int height, float zNear, float zFar)
{
	mat4f flip = mat4f::identity();
	flip(1, 1) = -1.f;
	flip(1, 2) = static_cast<float>(height);

	mat4f intrinsics = flip * intrinsic;
	const float L = 0.f;   // left
	const float R = static_cast<float>(width); // right
	const float B = 0.f;   // bottom
	const float T = static_cast<float>(height);// top

	// orthographic projection
	// https://msdn.microsoft.com/en-us/library/windows/desktop/bb205347(v=vs.85).aspx
	mat4f ortho = mat4f::zero();
	ortho(0, 0) = 2.0f / (R - L); ortho(0, 3) = (R + L) / (R - L);
	ortho(1, 1) = 2.0f / (T - B); ortho(1, 3) = (T + B) / (T - B);
	ortho(2, 2) = 1.0f / (zFar - zNear); ortho(2, 3) = zNear / (zNear - zFar);
	ortho(3, 3) = 1.f;

	// perspective projection
	mat4f tp = mat4f::zero();
	tp(0, 0) = intrinsics(0, 0); tp(0, 1) = intrinsics(0, 1); tp(0, 2) = -(width - 1 - intrinsics(0, 2));
	tp(1, 1) = intrinsics(1, 1); tp(1, 2) = -(height - 1 - intrinsics(1, 2));
	tp(2, 2) = zNear + zFar; tp(2, 3) = -zNear * zFar;
	tp(3, 2) = 1.0f;

	mat4f proj = ortho * tp;

	return proj;
}
