
#include "stdafx.h"

#include "scannedScene.h"
#include "imageHelper.h"
#include "normals.h"
#include "globalAppState.h"

#include "mLibFLANN.h"
#include "omp.h"


void ScannedScene::findKeyPoints()
{
	m_keyPoints.clear();
	const float depthSigmaD = GAS::get().s_depthFilterSigmaD;
	const float depthSigmaR = GAS::get().s_depthFilterSigmaR;
	const unsigned int maxNumKeyPoints = GAS::get().s_maxNumKeysPerFrame;
	const float minResponse = GAS::get().s_responseThresh;

	for (size_t sensorIdx = 0; sensorIdx < m_sds.size(); sensorIdx++) {
		SensorData* sd = m_sds[sensorIdx];
		const mat4f intrinsicInv = sd->m_calibrationDepth.m_intrinsic.getInverse();

		for (size_t i = 0; i < m_frameIndices[sensorIdx].size(); i++) {
			unsigned int imageIdx = m_frameIndices[sensorIdx][i];
		//const std::vector<unsigned int> imageIndices = { 34, 1276 }; //debugging
		//for (unsigned int imageIdx : imageIndices) { //debugging
			ColorImageR8G8B8 c = sd->computeColorImage(imageIdx);
			DepthImage32 d = sd->computeDepthImage(imageIdx);
			DepthImage32 dErode = d; ImageHelper::erode(dErode, 4);
			DepthImage32 dfilt = d; ImageHelper::bilateralFilter(dfilt, depthSigmaD, depthSigmaR);
			PointImage normalImage = SensorData::computeNormals(intrinsicInv, dfilt);

			//{ //debugging
			//	ColorImageR32G32B32 visNormalImage = normalImage;
			//	for (auto& p : visNormalImage) p.value = (p.value + 1.0f) / 0.5f; //scale [-1,1] to [0,1]
			//	FreeImageWrapper::saveImage("_depth.png", ColorImageR32G32B32(dfilt));
			//	FreeImageWrapper::saveImage("_normal.png", visNormalImage);
			//	const mat4f& camToWorld = sd->m_frames[imageIdx].getCameraToWorld();
			//	PointCloudf pc, pcw;
			//	for (const auto& p : d) {
			//		const vec3f& n = normalImage(p.x, p.y);
			//		if (n.x != -std::numeric_limits<float>::infinity()) {
			//			const vec3f c = intrinsicInv*vec3f(p.x*p.value, p.y*p.value, p.value);
			//			pc.m_points.push_back(c);
			//			pc.m_normals.push_back(n);
			//			pcw.m_points.push_back(camToWorld * c);
			//			pcw.m_normals.push_back((camToWorld.getRotation() * n).getNormalized());
			//		}
			//	}
			//	PointCloudIOf::saveToFile("_pcCam.ply", pc);
			//	PointCloudIOf::saveToFile("_pcWorld.ply", pcw);
			//	std::cout << "waiting..." << std::endl;
			//	getchar();
			//} //debugging

			//float sigmaD = 2.0f;	
			//float sigmaR = 0.1f;
			//sigmaD = 10.0f;
			//sigmaR = 10.0f;
			//FreeImageWrapper::saveImage("_before.png", c);
			//ImageHelper::bilateralFilter(c, sigmaD, sigmaR);
			//FreeImageWrapper::saveImage("_after.png", c);
			//std::cout << "here" << std::endl;
			//getchar();

			const mat4f& camToWorld = sd->m_frames[imageIdx].getCameraToWorld();

			std::vector<KeyPoint> rawKeyPoints = KeyPointFinder::findKeyPoints(vec2ui(sensorIdx, imageIdx), c, maxNumKeyPoints, minResponse);

			//MeshDataf md;
			size_t validKeyPoints = 0;
			for (KeyPoint& rawkp : rawKeyPoints) {
				const unsigned int padding = 50;	//don't take keypoints in the padding region of the image
				const vec2ui loc = math::round(rawkp.m_pixelPos);
				const vec2i dloc = math::round(vec2f((float)loc.x*(sd->m_depthWidth - 1) / (float)(sd->m_colorWidth - 1),
					(float)loc.y*(sd->m_depthHeight - 1) / (float)(sd->m_colorHeight-1)));
				const unsigned int dpad = std::round((float)padding*(sd->m_depthWidth - 1) / (float)(sd->m_colorWidth - 1));
				if (dErode.isValid(dloc) && dErode.isValidCoordinate(dloc + dpad) && dErode.isValidCoordinate(dloc - dpad)) {
					KeyPoint kp = rawkp;
					kp.m_depth = d(dloc);
					kp.m_imageIdx = (unsigned int)imageIdx;
					kp.m_sensorIdx = (unsigned int)sensorIdx;
					////kp.m_pixelPos = vec2f(rawkp.x, rawkp.y);
					kp.m_pixelPos = vec2f(loc);
					//kp.m_size = rawkp.m_size;
					//kp.m_angle = rawkp.m_angle;
					//kp.m_octave = rawkp.m_octave;
					//kp.m_scale = rawkp.m_scale;
					//kp.m_response = rawkp.m_response;
					//kp.m_opencvPackOctave = rawkp.m_opencvPackOctave;

					const vec3f cameraPos = (intrinsicInv*vec4f(kp.m_pixelPos.x*kp.m_depth, kp.m_pixelPos.y*kp.m_depth, kp.m_depth, 0.0f)).getVec3();
					kp.m_worldPos = camToWorld * cameraPos;

					const vec3f normal = normalImage(dloc); 
					if (normal.x == -std::numeric_limits<float>::infinity()) continue;
					kp.m_worldNormal = (camToWorld.getRotation() * normal).getNormalized();

					validKeyPoints++;

					m_keyPoints.push_back(kp);

					//if (imageIdx == 0) md.merge(Shapesf::sphere(0.01f, vec3f(kp.m_worldPos), 10, 10, vec4f(1.0f, 0.0f, 0.0f, 1.0f)).computeMeshData());
				}
			}

			std::cout << "\r" << "image: " << sensorIdx << "|" << imageIdx  << " found " << validKeyPoints << " keypoints";
			//if (imageIdx == 0) MeshIOf::saveToFile("test.ply", md);
			//if (imageIdx == 50) break;

			if (GAS::get().s_maxNumImages > 0 && imageIdx + 1 >= GAS::get().s_maxNumImages) break;
		}
	}
	std::cout << std::endl;
	std::cout << "total " << m_keyPoints.size() << " key points" << std::endl;
}

void ScannedScene::negativeKeyPoints()
{
	m_keyPointNegatives.clear();

	for (size_t i = 0; i < m_keyPointMatches.size(); i++) {
		KeyPointMatch noMatch;
		noMatch.m_kp0 = m_keyPointMatches[i].m_kp0;		//first key point is the same

		//search for non-matching keypoint in a different image.
		size_t tries = 0;
		for (;; tries++) {
			if (tries > 1000) throw MLIB_EXCEPTION("too many tries to find a negative match...");

			unsigned int idx = math::randomUniform(0u, (unsigned int)m_keyPoints.size() - 1);
			if (idx >= (unsigned int)m_keyPoints.size()) {
				std::cout << __FUNCTION__ << " ERROR: " << idx << "\t" << m_keyPoints.size() << std::endl;
			}
			noMatch.m_kp1 = m_keyPoints[idx];
			float dist = (noMatch.m_kp0.m_worldPos - noMatch.m_kp1.m_worldPos).length();
			if (noMatch.m_kp0.isSameImage(noMatch.m_kp1) || dist <= GAS::get().s_matchThresh) continue;	//invalid; either same image or within the threshold
			else {
				//found a good match
				noMatch.m_offset = vec2f(0.0f);
				m_keyPointNegatives.push_back(noMatch);
				break;
			}
		}
	}

	if (m_keyPointMatches.size() != m_keyPointNegatives.size()) throw MLIB_EXCEPTION("something went wrong here... positive and negative match counts should be the same");
}

void ScannedScene::matchKeyPoints()
{
	const float radius = GAS::get().s_matchThresh;
	const unsigned int maxK = 5;
	const unsigned int maxNumMatchesPerScene = GAS::get().s_maxNumMatchesPerScene;

	unsigned int currKeyPoint = 0;
	std::vector<std::vector<NearestNeighborSearchFLANNf*>> nns(m_sds.size());
	std::vector<std::vector<unsigned int>> nn_offsets(m_sds.size());
	for (size_t sensorIdx = 0; sensorIdx < nns.size(); sensorIdx++) {
		nns[sensorIdx].resize(m_sds[sensorIdx]->m_frames.size(), nullptr);
		nn_offsets[sensorIdx].resize(m_sds[sensorIdx]->m_frames.size(), 0);
		for (size_t i = 0; i < m_frameIndices[sensorIdx].size(); i++) {
			const unsigned int frameIdx = m_frameIndices[sensorIdx][i];
			nn_offsets[sensorIdx][frameIdx] = currKeyPoint;

			std::vector<float> points;
			for (;	currKeyPoint < (UINT)m_keyPoints.size() &&
				m_keyPoints[currKeyPoint].m_sensorIdx == sensorIdx &&
				m_keyPoints[currKeyPoint].m_imageIdx == frameIdx; currKeyPoint++) {
				points.push_back(m_keyPoints[currKeyPoint].m_worldPos.x);
				points.push_back(m_keyPoints[currKeyPoint].m_worldPos.y);
				points.push_back(m_keyPoints[currKeyPoint].m_worldPos.z);
			}

			if (points.size())	{
				nns[sensorIdx][frameIdx] = new NearestNeighborSearchFLANNf(4 * maxK, 1);
				nns[sensorIdx][frameIdx]->init(points.data(), (unsigned int)points.size() / 3, 3, maxK);
			}
		}
	}


	for (size_t keyPointIdx = 0; keyPointIdx < m_keyPoints.size(); keyPointIdx++) {

		if (keyPointIdx % 100 == 0)
			std::cout << "\rmatching keypoint " << keyPointIdx << " out of " << m_keyPoints.size();

		KeyPoint& kp = m_keyPoints[keyPointIdx];
		const float* query = (const float*)&kp.m_worldPos;

		const size_t sensorIdx = kp.m_sensorIdx;
		const size_t frameIdx = kp.m_imageIdx;

		for (size_t sensorIdx_dst = sensorIdx; sensorIdx_dst < nns.size(); sensorIdx_dst++) {
			if (maxNumMatchesPerScene > 0 && m_keyPointMatches.size() > maxNumMatchesPerScene) break;

			size_t frameIdx_dst = 0;
			if (sensorIdx_dst == sensorIdx) frameIdx_dst = frameIdx + 1;

			for (; frameIdx_dst < nns[sensorIdx].size(); frameIdx_dst++) {
				if (maxNumMatchesPerScene > 0 && m_keyPointMatches.size() > maxNumMatchesPerScene) break;

				auto* nn = nns[sensorIdx_dst][frameIdx_dst];
				if (nn == nullptr) continue;

				size_t numMatches = 0;
				//std::vector<unsigned int> res = nn->fixedRadius(query, maxK, radius);
				auto resPair = nn->fixedRadiusDist(query, maxK, radius);
				//auto resDist = nn->getDistances((UINT)res.size());
				std::sort(resPair.begin(), resPair.end(), [](const std::pair<unsigned int, float> &left, const std::pair<unsigned int, float> &right) {
					return fabs(left.second) < fabs(right.second);
				});
				for (size_t j = 0; j < resPair.size(); j++) {
					const KeyPoint& kp0 = m_keyPoints[keyPointIdx];
					const KeyPoint& kp1 = m_keyPoints[resPair[j].first + nn_offsets[sensorIdx_dst][frameIdx_dst]];

					//check world normals
					float angleDist = std::acos(kp0.m_worldNormal | kp1.m_worldNormal);
					if (angleDist > 1.75f) continue; //ignore with world normals > ~100 degrees apart
					
					////debugging
					//const float degrees = math::radiansToDegrees(angleDist);
					//if (degrees > 90) { 
					//	KeyPointMatch m;
					//	m.m_kp0 = kp0;
					//	m.m_kp1 = kp1;
					//	MatchVisualization mv;
					//	mv.visulizeMatches(m_sds, std::vector<KeyPointMatch>(1, m), 10);
					//	m_sds[kp0.m_sensorIdx]->saveToPointCloud("frame0.ply", kp0.m_imageIdx);
					//	MeshDataf keymesh = Shapesf::sphere(0.03f, kp0.m_worldPos, 10, 10, vec4f(1.0f, 0.0f, 0.0f, 1.0f)).computeMeshData();
					//	keymesh.merge(Shapesf::cylinder(kp0.m_worldPos, kp0.m_worldPos + kp0.m_worldNormal * 0.2f, .02f, 10, 10, vec4f(1.0f, 0.0f, 0.0f, 1.0f)).computeMeshData());
					//	MeshIOf::saveToFile("frame0_key.ply", keymesh);
					//	m_sds[kp1.m_sensorIdx]->saveToPointCloud("frame1.ply", kp1.m_imageIdx);
					//	keymesh = Shapesf::sphere(0.03f, kp1.m_worldPos, 10, 10, vec4f(0.0f, 1.0f, 0.0f, 1.0f)).computeMeshData();
					//	keymesh.merge(Shapesf::cylinder(kp1.m_worldPos, kp1.m_worldPos + kp1.m_worldNormal * 0.2f, .02f, 10, 10, vec4f(0.0f, 1.0f, 0.0f, 1.0f)).computeMeshData());
					//	MeshIOf::saveToFile("frame1_key.ply", keymesh);
					//	int a = 5;
					//} //debugging

					KeyPointMatch m;
					m.m_kp0 = kp0;
					m.m_kp1 = kp1;

					mat4f worldToCam_kp1 = m_sds[m.m_kp1.m_sensorIdx]->m_frames[m.m_kp1.m_imageIdx].getCameraToWorld().getInverse();
					vec3f p = worldToCam_kp1 * m.m_kp0.m_worldPos;
					p = m_sds[m.m_kp1.m_sensorIdx]->m_calibrationDepth.cameraToProj(p);
					m.m_offset = m.m_kp1.m_pixelPos - vec2f(p.x, p.y);
					m_keyPointMatches.push_back(m);
					numMatches++;
					break;

					//std::cout << "orig: " << m.m_kp1.m_pixelPos << std::endl;
					//std::cout << "repr: " << vec2f(p.x, p.y) << std::endl;
					//std::cout << m.m_offset << std::endl;

					//std::cout << "match between: " << std::endl;
					//std::cout << m.m_kp0;
					//std::cout << m.m_kp1;

					//std::cout << "dist " << sensorIdx << ":\t" << (m.m_kp0.m_worldPos - m.m_kp1.m_worldPos).length() << std::endl;
					//std::cout << std::endl;
					//int a = 5;
				}
			}
		}
	}

	std::cout << std::endl;


	//clean up our mess...
	for (size_t sensorIdx = 0; sensorIdx < nns.size(); sensorIdx++) {
		for (size_t frameIdx = 0; frameIdx < nns[sensorIdx].size(); frameIdx++) {
			SAFE_DELETE(nns[sensorIdx][frameIdx]);
		}
	}


	std::cout << "Total matches found: " << m_keyPointMatches.size() << std::endl;
}

void ScannedScene::saveImages(const std::string& outPath) const
{
	if (!util::directoryExists(outPath)) {
		util::makeDirectory(outPath);
	}

	const bool bSaveNormals = !m_normals.empty();

	const unsigned int outWidth = GAS::get().s_outWidth;
	const unsigned int outHeight = GAS::get().s_outHeight;
	const unsigned int maxNumSensors = std::min((unsigned int)m_sds.size(), GAS::get().s_maxNumSensFiles);
	for (unsigned int sensorIdx = 0; sensorIdx < maxNumSensors; sensorIdx++) {
		SensorData* sd = m_sds[sensorIdx];

		for (size_t i = 0; i < m_frameIndices[sensorIdx].size(); i++) {
			const unsigned int imageIdx = m_frameIndices[sensorIdx][i];
			ColorImageR8G8B8 c = sd->computeColorImage(imageIdx);
			c.resize(outWidth, outHeight);
			//DepthImage32 d = sd->computeDepthImage(imageIdx);
			//TODO depth image?

			char s[256];
			sprintf(s, "color-%02d-%06d.jpg", (UINT)sensorIdx, (UINT)imageIdx);
			const std::string outFileColor = outPath + "/" + std::string(s);

			std::cout << "\r" << outFileColor;
			FreeImageWrapper::saveImage(outFileColor, c);

			if (GAS::get().s_maxNumImages > 0 && imageIdx + 1 >= GAS::get().s_maxNumImages) break;
		}
		std::cout << std::endl;
	}
}

void ScannedScene::computeNormals(NORMAL_TYPE type)
{
	switch (type) 
	{
	case DEPTH_NORMALS:
		computeDepthNormals();
		break;
	case MESH_NORMALS:
		computeMeshNormals();
		break;
	default:
		//should never come here
		break;
	}
}

void ScannedScene::computeDepthNormals()
{
	std::cout << "computing depth normals... ";
	const unsigned int maxNumSensors = std::min((unsigned int)m_sds.size(), GAS::get().s_maxNumSensFiles);
	m_normals.resize(maxNumSensors);
	NormalExtractor extractor(GAS::get().s_outWidth, GAS::get().s_outHeight);
	for (unsigned int sensorIdx = 0; sensorIdx < maxNumSensors; sensorIdx++) {
		extractor.computeDepthNormals(m_sds[sensorIdx], GAS::get().s_depthFilterSigmaD, GAS::get().s_depthFilterSigmaR, m_normals[sensorIdx]);
	}
	std::cout << "done!" << std::endl;
}

void ScannedScene::computeMeshNormals()
{
	std::cout << "computing mesh normals... ";
	//load mesh
	const std::string meshFile = m_path + "/" + m_name + "_vh_clean.ply"; //highest-res
	TriMeshf mesh = TriMeshf(MeshIOf::loadFromFile(meshFile));
	mesh.computeNormals();
	const unsigned int maxNumSensors = std::min((unsigned int)m_sds.size(), GAS::get().s_maxNumSensFiles);
	m_normals.resize(maxNumSensors);
	std::vector<NormalExtractor*> extractors(omp_get_max_threads());
	for (unsigned int i = 0; i < extractors.size(); i++) extractors[i] = new NormalExtractor(GAS::get().s_outWidth, GAS::get().s_outHeight);
#pragma omp parallel for
	for (int sensorIdx = 0; sensorIdx < (int)maxNumSensors; sensorIdx++) {
		const int thread = omp_get_thread_num();
		NormalExtractor& extractor = *extractors[thread];
		extractor.computeMeshNormals(m_sds[sensorIdx], mesh, GAS::get().s_renderDepthMin, GAS::get().s_renderDepthMax, m_normals[sensorIdx]);
	}

	for (unsigned int i = 0; i < extractors.size(); i++)
		SAFE_DELETE(extractors[i]);
	std::cout << "done!" << std::endl;
}

template<>
struct std::hash<vec2ui> : public std::unary_function < vec2ui, size_t > {
	size_t operator()(const vec2ui& v) const {
		//TODO larger prime number (64 bit) to match size_t
		const size_t p0 = 73856093;
		const size_t p1 = 19349669;
		//const size_t p2 = 83492791;
		const size_t res = ((size_t)v.x * p0) ^ ((size_t)v.y * p1);// ^ ((size_t)v.z * p2);
		return res;
	}
};

void ScannedScene::saveNormalImages(const std::string& outPath) const
{
	if (m_normals.empty()) {
		std::cout << "no normals to save" << std::endl;
		return;
	}

	if (!util::directoryExists(outPath)) {
		util::makeDirectory(outPath);
	}
	std::string rawSrcPath = GAS::get().s_srcPath + "/" + m_name + "/data/";
	if (!util::directoryExists(rawSrcPath)) rawSrcPath = "//falas/Matterport/v1/" + m_name + "/matterport_depth_images/";
	if (!util::directoryExists(rawSrcPath)) throw MLIB_EXCEPTION("raw src path (" + rawSrcPath + ") does not exist");
	Directory dir(rawSrcPath);
	std::vector<std::string> baseFiles = dir.getFilesWithSuffix("_d0_0.png");
	for (auto& f : baseFiles) f = util::replace(f, "_d0_0.png", "");

	std::unordered_map<vec2ui, std::string> nameMap;
	{
		std::vector<unsigned int> imageIdxs(3);
		for (size_t fidx = 0; fidx < baseFiles.size(); fidx++) {
			const std::string& f = baseFiles[fidx];
			for (unsigned int camIdx = 0; camIdx < 3; camIdx++) {
				for (unsigned int poseIdx = 0; poseIdx < 6; poseIdx++) {
					nameMap[vec2ui(camIdx, imageIdxs[camIdx])] = f + "_n" + std::to_string(camIdx) + "_" + std::to_string(poseIdx);
					MLIB_ASSERT(util::fileExists(rawSrcPath + f + "_d" + std::to_string(camIdx) + "_" + std::to_string(poseIdx) + ".png"));
					imageIdxs[camIdx]++;
				}
			}
		}
	}

	const unsigned int maxNumSensors = std::min(GAS::get().s_maxNumSensFiles, (unsigned int)m_sds.size());
	for (unsigned int sensorIdx = 0; sensorIdx < maxNumSensors; sensorIdx++) {
		SensorData* sd = m_sds[sensorIdx];
		const mat4f depthIntrinsicInv = sd->m_calibrationDepth.m_intrinsic.getInverse();

		for (size_t imageIdx = 0; imageIdx < sd->m_frames.size(); imageIdx++) {

			//char s[256];
			//sprintf(s, "normal-%02d-%06d.png", (UINT)sensorIdx, (UINT)imageIdx);
			//const std::string outFileNormal = outPath + "/" + std::string(s);
			const auto it = nameMap.find(vec2ui(sensorIdx, (unsigned int)imageIdx));
			if (it == nameMap.end()) throw MLIB_EXCEPTION("no name for sens " + std::to_string(sensorIdx) + ", image " + std::to_string(imageIdx));
			const std::string outFileNormal = outPath + "/" + it->second + ".png";

			std::cout << "\r" << outFileNormal;
			NormalExtractor::saveNormalImage(outFileNormal, m_normals[sensorIdx][imageIdx]);

			if (GAS::get().s_maxNumImages > 0 && imageIdx + 1 >= GAS::get().s_maxNumImages) break;
		}
		std::cout << std::endl;
	}
}

void ScannedScene::debug()
{
	const std::string file = "W:/data/matterport/derived_data/17DRP5sb8fy/normals_depth/00ebbf3782c64d74aaf7dd39cd561175_n1_2.png";
	PointImage normal;
	NormalExtractor::loadNormalImage(file, normal);

	ColorImageR8G8B8 color(normal.getDimensions());
	for (unsigned int y = 0; y < normal.getHeight(); y++) {
		for (unsigned int x = 0; x < normal.getWidth(); x++) {
			vec3f n = normal(x, y);
			if (n.x == -std::numeric_limits<float>::infinity()) color(x, y) = vec3uc(0, 0, 0);
			else {
				n = (n + 1.0f) * 0.5f * 255.0f;
				color(x, y) = vec3uc(n);
			}
		}
	}
	FreeImageWrapper::saveImage("debug.png", color);
	int a = 5;
}

void ScannedScene::debugMatch() const
{
	const vec2ui im0(0, 8);
	const vec2ui im1(1, 234);
	const vec2ui loc0(801, 896);
	const vec2ui loc1(546, 457);

	const vec2ui imn(2, 141);
	const vec2ui locn(994, 723);

	//12200
	//const vec2ui im0(0, 124);
	//const vec2ui im1(1, 507);
	//const vec2ui loc0(188, 474);
	//const vec2ui loc1(416, 470);

	DepthImage32 d0 = m_sds[im0.x]->computeDepthImage(im0.y);
	DepthImage32 d1 = m_sds[im1.x]->computeDepthImage(im1.y);
	KeyPoint kp0;
	kp0.m_sensorIdx = im0.x; kp0.m_imageIdx = im0.y;
	kp0.m_pixelPos = vec2f(loc0); kp0.m_depth = d0(loc0);
	kp0.m_worldPos = m_sds[im0.x]->m_frames[im0.y].getCameraToWorld() * (m_sds[im0.x]->m_calibrationDepth.m_intrinsic.getInverse() *vec3f(kp0.m_pixelPos.x*kp0.m_depth, kp0.m_pixelPos.y*kp0.m_depth, kp0.m_depth));
	KeyPoint kp1;
	kp1.m_sensorIdx = im1.x; kp1.m_imageIdx = im1.y;
	kp1.m_pixelPos = vec2f(loc1); kp1.m_depth = d1(loc1);
	kp1.m_worldPos = m_sds[im1.x]->m_frames[im1.y].getCameraToWorld() * (m_sds[im1.x]->m_calibrationDepth.m_intrinsic.getInverse() *vec3f(kp1.m_pixelPos.x*kp1.m_depth, kp1.m_pixelPos.y*kp1.m_depth, kp1.m_depth));

	//compute world normals
	const float depthSigmaD = GAS::get().s_depthFilterSigmaD;
	const float depthSigmaR = GAS::get().s_depthFilterSigmaR;
	ImageHelper::bilateralFilter(d0, depthSigmaD, depthSigmaR);
	ImageHelper::bilateralFilter(d1, depthSigmaD, depthSigmaR);
	PointImage normal0 = SensorData::computeNormals(m_sds[im0.x]->m_calibrationDepth.m_intrinsic.getInverse(), d0);
	PointImage normal1 = SensorData::computeNormals(m_sds[im1.x]->m_calibrationDepth.m_intrinsic.getInverse(), d1);
	const vec3f camNormal0 = normal0(loc0);		const vec3f worldNormal0 = m_sds[im0.x]->m_frames[im0.y].getCameraToWorld().getRotation() * camNormal0;
	const vec3f camNormal1 = normal1(loc1);		const vec3f worldNormal1 = m_sds[im1.x]->m_frames[im1.y].getCameraToWorld().getRotation() * camNormal1;
	const float angleDist = std::acos(worldNormal0 | worldNormal1);

	KeyPointMatch m;
	m.m_kp0 = kp0;
	m.m_kp1 = kp1;
	MatchVisualization mv;
	mv.visulizeMatches(m_sds, std::vector<KeyPointMatch>(1, m), 10, 1, true);
	mv.visulizeMatches3D(m_sds, std::vector<KeyPointMatch>(1, m), 10);

	ColorImageR8G8B8 ca = m_sds[im0.x]->computeColorImage(im0.y); ca.resize(640, 512);
	ColorImageR8G8B8 cp = m_sds[im1.x]->computeColorImage(im1.y); cp.resize(640, 512);
	ColorImageR8G8B8 cn = m_sds[imn.x]->computeColorImage(imn.y); cn.resize(640, 512);
	const unsigned int dim = 65; const int radius = 32;
	ColorImageR8G8B8 patchAnc(dim, dim);
	ColorImageR8G8B8 patchPos(dim, dim);
	ColorImageR8G8B8 patchNeg(dim, dim);
	for (int y = -radius; y <= radius; y++) {
		for (int x = -radius; x <= radius; x++) {
			patchAnc(x + radius, y + radius) = ca(math::round(vec2f(loc0)*0.5f) + vec2i(x, y));
			patchPos(x + radius, y + radius) = cp(math::round(vec2f(loc1)*0.5f) + vec2i(x, y));
			patchNeg(x + radius, y + radius) = cn(math::round(vec2f(locn)*0.5f) + vec2i(x, y));
		}
	}
	FreeImageWrapper::saveImage("patch_anc.png", patchAnc);
	FreeImageWrapper::saveImage("patch_pos.png", patchPos);
	FreeImageWrapper::saveImage("patch_neg.png", patchNeg);

	//DepthImage32 de = d1; ImageHelper::erode(de, 2);
	//FreeImageWrapper::saveImage("_depth.png", ColorImageR32G32B32(d1));
	//FreeImageWrapper::saveImage("_depth-erode2.png", ColorImageR32G32B32(de));
	//de = d1; ImageHelper::erode(de, 4);
	//FreeImageWrapper::saveImage("_depth-erode4.png", ColorImageR32G32B32(de));

	std::cout << "waiting..." << std::endl; getchar();
}
