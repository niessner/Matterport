
#include "stdafx.h"

#include "scannedScene.h"
#include "imageHelper.h"

//#define DEBUG_IMAGES_POTENTIAL
//#define DEBUG_IMAGES_REPROJ

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

ReprojError ScannedScene::computeReprojection(size_t maxNumFramePairSamples, size_t maxNumSampleTries, float maxDepth, std::vector<float>& sampledCameraDists)
{
	size_t numFramesPerSens = m_sds.front()->m_frames.size();
	size_t numFrames = numFramesPerSens * m_sds.size();
	const size_t numFramePairs = numFrames*(numFrames - 1) / 2;
	if (maxNumFramePairSamples > numFramePairs) { //should not happen
		std::cout << "warning: not enough frames (" << numFrames << ") for " << maxNumFramePairSamples << " pair samples" << std::endl;
		maxNumFramePairSamples = numFramePairs;
	}
	const unsigned int width = m_sds.front()->m_depthWidth;
	const unsigned int height = m_sds.front()->m_depthHeight;

#ifdef FIND_FAR_CAMERAS //doesn't actually enforce a distribution, just looks for far away cameras..
	MLIB_ASSERT(m_sds.size() == 1); //only for scannet
	const float avgCameraDist = 2.62f; const float stdCameraDist = 1.52f;
	std::vector< std::pair<vec2ui, float> > framePairCameraDists;
	for (unsigned int i = 0; i < numFrames; i++) {
		const mat4f& t0 = m_sds[0]->m_frames[i].getCameraToWorld();
		if (t0[0] == -std::numeric_limits<float>::infinity()) continue;
		for (unsigned int j = i + 1; j < numFrames; j++) {
			const mat4f& t1 = m_sds[0]->m_frames[j].getCameraToWorld();
			if (t1[0] == -std::numeric_limits<float>::infinity()) continue;
			const float dist =  vec3f::dist(t0.getTranslation(), t1.getTranslation());
			if (dist < 6.0f)
				framePairCameraDists.push_back(std::make_pair(vec2ui(i, j), dist));
		}
	}
	std::sort(framePairCameraDists.begin(), framePairCameraDists.end(), [](const std::pair<vec2ui, float> &left, const std::pair<vec2ui, float> &right) {
		return fabs(left.second) > fabs(right.second); //sort largest camera dist first
	});
	unsigned int curFramePairIdx = 0;
	unsigned int skip = 1;
	if (numFrames > 2000 && numFrames < 3000) skip = 10;
	else if (numFrames > 3000) skip = 20;
#endif

	ReprojError err;
	std::unordered_set<vec2ui> framePairs; 
	size_t numTries = 0;
	unsigned int numFound = 0;
	while (numFound < maxNumFramePairSamples && numTries < maxNumSampleTries) {
		numTries++;
#ifdef FIND_FAR_CAMERAS
		if (numFound == 0 && numTries % 100 == 0) std::cout << "\rtry " << numTries << "\t" << curFramePairIdx << "\t" << framePairCameraDists[curFramePairIdx].second;
		const unsigned int f0 = framePairCameraDists[curFramePairIdx].first.x;
		const unsigned int f1 = framePairCameraDists[curFramePairIdx].first.y;
		curFramePairIdx+=skip;
		if (curFramePairIdx >= framePairCameraDists.size()) {
			std::cout << "exhausted all frame pair cameras (" << framePairCameraDists.size() << ") in " << numTries << " tries" << std::endl;
			break;
		}
		//if (maxNumFramePairSamples - numFound < framePairCameraDists.size()-curFramePairIdx) { //can still reject while maintaining requested #samples
		//	//if (!sampledCameraDists.empty()) {
		//	//	const float curAvgCameraDist = std::accumulate(sampledCameraDists.begin(), sampledCameraDists.end())
		//	//}
		//}
#else
		const unsigned int f0 = math::randomUniform(0u, (unsigned int)numFrames - 1);
		const unsigned int f1 = math::randomUniform(0u, (unsigned int)numFrames - 1);
#endif
		const vec2ui frames = f0 < f1 ? vec2ui(f0, f1) : vec2ui(f1, f0);
		if (f0 == f1 || framePairs.find(frames) != framePairs.end()) continue; //not a valid sample
		framePairs.insert(frames);

		vec2ui sensFrameIdx0(0, f0), sensFrameIdx1(0, f1);
		if (m_sds.size() > 1) {
			sensFrameIdx0.x = f0 / (unsigned int)numFramesPerSens;
			sensFrameIdx0.y = f0 % (unsigned int)numFramesPerSens;
			sensFrameIdx1.x = f1 / (unsigned int)numFramesPerSens;
			sensFrameIdx1.y = f1 % (unsigned int)numFramesPerSens;
		}
		const mat4f& t0 = m_sds[sensFrameIdx0.x]->m_frames[sensFrameIdx0.y].getCameraToWorld();
		const mat4f& t1 = m_sds[sensFrameIdx1.x]->m_frames[sensFrameIdx1.y].getCameraToWorld();
		if (t0[0] == -std::numeric_limits<float>::infinity() || t1[0] == -std::numeric_limits<float>::infinity()) 
			continue;
		const mat4f transform0to1 = t1.getInverse() * t0;
		const mat4f depthIntrinsics0 = m_sds[sensFrameIdx0.x]->m_calibrationDepth.m_intrinsic;
		const mat4f depthIntrinsicsInv0 = m_sds[sensFrameIdx0.x]->m_calibrationDepth.m_intrinsic.getInverse();
		const mat4f depthIntrinsics1 = m_sds[sensFrameIdx1.x]->m_calibrationDepth.m_intrinsic;
		const mat4f depthIntrinsicsInv1 = m_sds[sensFrameIdx1.x]->m_calibrationDepth.m_intrinsic.getInverse();
		//see if there is any potential overlap between the frames
		if (!Reprojection::hasCameraFrustumIntersection(t0, depthIntrinsics0, t1, depthIntrinsics1, width, height, maxDepth))
			continue; //no intersection

		DepthImage32 depth0 = m_sds[sensFrameIdx0.x]->computeDepthImage(sensFrameIdx0.y);
		DepthImage32 depth1 = m_sds[sensFrameIdx1.x]->computeDepthImage(sensFrameIdx1.y);

		if (!Reprojection::checkValidReprojection(depth0, depth1, depthIntrinsics0, depthIntrinsics1, depthIntrinsicsInv0, depthIntrinsicsInv1, transform0to1, maxDepth))
			continue; //no valid projections

		ColorImageR32 intensity0, intensity1;
		{
			ColorImageR8G8B8 color0 = m_sds[sensFrameIdx0.x]->computeColorImage(sensFrameIdx0.y);
			ColorImageR8G8B8 color1 = m_sds[sensFrameIdx1.x]->computeColorImage(sensFrameIdx1.y);
			//MLIB_ASSERT((float)color0.getWidth() / depth0.getWidth() == (float)color0.getHeight() / depth0.getHeight());
			if (color0.getWidth() != depth0.getWidth()) { //assumes extrinsics == identity
				color0.resize(depth0.getWidth(), depth0.getHeight());
				color1.resize(depth1.getWidth(), depth1.getHeight());
			}
			intensity0 = ImageHelper::convertToGrayscale(color0);
			intensity1 = ImageHelper::convertToGrayscale(color1);
		}

#ifdef DEBUG_IMAGES_POTENTIAL
		FreeImageWrapper::saveImage("depth0.png", ColorImageR32G32B32(depth0));
		FreeImageWrapper::saveImage("depth1.png", ColorImageR32G32B32(depth1));
		FreeImageWrapper::saveImage("intensity0.png", intensity0);
		FreeImageWrapper::saveImage("intensity1.png", intensity1);
		{
			PointCloudf pp0, pp1;
			for (const auto& p : depth0) {
				const vec3f p0 = Reprojection::computeCameraSpacePosition(depth0, depthIntrinsicsInv0, p.x, p.y);
				const vec3f p1 = Reprojection::computeCameraSpacePosition(depth1, depthIntrinsicsInv1, p.x, p.y);
				if (p0.x != -std::numeric_limits<float>::infinity()) pp0.m_points.push_back(t0 * p0);
				if (p1.x != -std::numeric_limits<float>::infinity()) pp1.m_points.push_back(t1 * p1);
			}
			PointCloudIOf::saveToFile("pc0.ply", pp0);
			PointCloudIOf::saveToFile("pc1.ply", pp1);
		}
		{
			std::vector<vec3f> frustumPoints; //everything in camera space
			frustumPoints.push_back(vec3f::origin);
			frustumPoints.push_back(depthIntrinsicsInv0 * (maxDepth * vec3f(0, 0, 1.0f)));
			frustumPoints.push_back(depthIntrinsicsInv0 * (maxDepth * vec3f(width - 1.0f, 0, 1.0f)));
			frustumPoints.push_back(depthIntrinsicsInv0 * (maxDepth * vec3f(width - 1.0f, height - 1.0f, 1.0f)));
			frustumPoints.push_back(depthIntrinsicsInv0 * (maxDepth * vec3f(0, height - 1.0f, 1.0f)));

			MeshDataf f0, f1;
			f0.merge(Shapesf::cylinder(t0 * frustumPoints[0], t0 * frustumPoints[1], 0.02f, 10, 10).computeMeshData());
			f0.merge(Shapesf::cylinder(t0 * frustumPoints[0], t0 * frustumPoints[2], 0.02f, 10, 10).computeMeshData());
			f0.merge(Shapesf::cylinder(t0 * frustumPoints[0], t0 * frustumPoints[3], 0.02f, 10, 10).computeMeshData());
			f0.merge(Shapesf::cylinder(t0 * frustumPoints[0], t0 * frustumPoints[4], 0.02f, 10, 10).computeMeshData());

			f1.merge(Shapesf::cylinder(t1 * frustumPoints[0], t1 * frustumPoints[1], 0.02f, 10, 10).computeMeshData());
			f1.merge(Shapesf::cylinder(t1 * frustumPoints[0], t1 * frustumPoints[2], 0.02f, 10, 10).computeMeshData());
			f1.merge(Shapesf::cylinder(t1 * frustumPoints[0], t1 * frustumPoints[3], 0.02f, 10, 10).computeMeshData());
			f1.merge(Shapesf::cylinder(t1 * frustumPoints[0], t1 * frustumPoints[4], 0.02f, 10, 10).computeMeshData());

			MeshIOf::saveToFile("f0.ply", f0);
			MeshIOf::saveToFile("f1.ply", f1);
		}
		std::cout << "found potential overlap" << std::endl;
		getchar();
#endif

		ImageHelper::bilateralFilter(depth0, 2.0f, 0.05f);
		ImageHelper::bilateralFilter(depth1, 2.0f, 0.05f);
		ImageHelper::gaussFilter(intensity0, 2.0f);
		ImageHelper::gaussFilter(intensity1, 2.0f);

		ReprojError e = Reprojection::computeReprojection(depth0, intensity0, depth1, intensity1, 
			depthIntrinsics0, depthIntrinsics1, depthIntrinsicsInv0, depthIntrinsicsInv1, transform0to1, maxDepth);
		if (e.numCorrs == 0) 
			continue; // no intersection
		////debugging
		//const bool first = !util::fileExists("cameras.csv");
		//std::ofstream ofs("cameras.csv", std::ios::app);
		//if (first) ofs << "attempt,camera distance,campos0,campos1" << std::endl;
		//ofs << numFound << "," << vec3f::dist(t0.getTranslation(), t1.getTranslation()) << "," << t0.getTranslation() << "," << t1.getTranslation() << std::endl;
		//ofs.close();
		////debugging
#ifdef FIND_FAR_CAMERAS
		sampledCameraDists.push_back(vec3f::dist(t0.getTranslation(), t1.getTranslation()));
		if (numFound == 0) {
			framePairCameraDists.erase(framePairCameraDists.begin(), framePairCameraDists.begin() + curFramePairIdx);
			unsigned int idx24 = 0, idx2 = 0, idx1 = 0, idx05 = 0;
			for (unsigned int i = 0; i < framePairCameraDists.size(); i++) {
				if (idx24 == 0 && framePairCameraDists[i].second < 2.4f) idx24 = i;
				if (idx2 == 0 && framePairCameraDists[i].second < 2.0f) idx2 = i;
				if (idx1 == 0 && framePairCameraDists[i].second < 1.0f) idx1 = i;
				if (idx05 == 0 && framePairCameraDists[i].second < 0.5f) idx05 = i;
			}
			std::cout << "\rnew set size = " << framePairCameraDists.size() << "\t(" << framePairCameraDists.front().second << "," << framePairCameraDists.back().second << ")" << std::endl;
			{
				std::vector< std::pair<vec2ui, float> > tmp; std::vector<unsigned int> indices;
				indices = std::vector<unsigned int>(std::min((unsigned int)maxNumFramePairSamples * 90, idx24)); //indices for dists > 2.4m
				for (unsigned int i = 0; i < indices.size(); i++) indices[i] = i;
				std::random_shuffle(indices.begin(), indices.end());
				for (const auto& idx : indices) tmp.push_back(framePairCameraDists[idx]);

				indices = std::vector<unsigned int>(std::min((unsigned int)maxNumFramePairSamples * 10, idx24 - idx2)); //indices for dists > 2m
				for (unsigned int i = idx24; i < idx24 + indices.size(); i++) indices[i - idx24] = i;
				std::random_shuffle(indices.begin(), indices.end());
				for (const auto& idx : indices) tmp.push_back(framePairCameraDists[idx]);

				indices = std::vector<unsigned int>(std::min((unsigned int)(maxNumFramePairSamples * 0.5f), idx1 - idx2)); //indices for dists > 1m
				std::random_shuffle(indices.begin(), indices.end());
				for (unsigned int i = idx2; i < idx2 + indices.size(); i++) indices[i - idx2] = i;
				for (const auto& idx : indices) tmp.push_back(framePairCameraDists[idx]);

				indices = std::vector<unsigned int>(std::min((unsigned int)(maxNumFramePairSamples * 0.1f), idx05 - idx1)); //indices for dists > 0.5m
				std::random_shuffle(indices.begin(), indices.end());
				for (unsigned int i = idx1; i < idx1 + indices.size(); i++) indices[i - idx1] = i;
				for (const auto& idx : indices) tmp.push_back(framePairCameraDists[idx]);

				framePairCameraDists = tmp;
			}
			std::cout << "trunc set size = " << framePairCameraDists.size() << std::endl;
			std::random_shuffle(framePairCameraDists.begin(), framePairCameraDists.end());
			curFramePairIdx = 0;
			skip = 1;
		}
#endif
		err += e;
		numFound++;
		std::cout << "\rfound " << numFound << " (" << numTries << " attempts)";
	}
	if (numTries == maxNumSampleTries) std::cout << "\rexhausted #tries, found " << numFound << " frame pairs" << std::endl;
	else std::cout << "\rfound " << numFound << " frame pairs in " << numTries << " tries" << std::endl;
	return err;
}

ReprojError Reprojection::computeReprojection(const DepthImage32& depth0, const ColorImageR32& intensity0, 
	const DepthImage32& depth1, const ColorImageR32& intensity1, 
	const mat4f& depthIntrinsics0, const mat4f& depthIntrinsics1, 
	const mat4f& depthIntrinsicsInv0, const mat4f& depthIntrinsicsInv1, const mat4f& transform0to1, float maxDepth)
{
	const unsigned int subsampleFactor = (depth0.getWidth() == 640) ? 8 : 16; //todo fix hack
	const float normalThresh = 0.95f;
	const float distThresh = 0.15f;
	const float colorThresh = 0.1f;
	const float depthMin = 0.4f;
	const float depthMax = maxDepth;
	const float INVALID = -std::numeric_limits<float>::infinity();

	ColorImageR32 gradmag0 = ImageHelper::computeGradientMagnitude(intensity0);
	ColorImageR32 gradmag1 = ImageHelper::computeGradientMagnitude(intensity1);

#ifdef DEBUG_IMAGES_REPROJ
	ImageHelper::computeImageStatistics(gradmag0);
	ImageHelper::computeImageStatistics(gradmag1);
	DepthImage32 debug(depth0.getWidth()/subsampleFactor, depth0.getHeight()/subsampleFactor);
	debug.setInvalidValue(INVALID); debug.setPixels(INVALID);
	PointCloudf pc0trans, pc1;
#endif
	ReprojError err;
#pragma omp parallel for
	for (int _y = 0; _y < (int)depth0.getHeight(); _y+=subsampleFactor) {
		unsigned int y = (unsigned int)_y;
		for (unsigned int x = 0; x < depth0.getWidth(); x+=subsampleFactor) {
			const float d0 = depth0(x, y);
			const vec3f p0 = computeCameraSpacePosition(depth0, depthIntrinsicsInv0, x, y);
			const vec3f n0 = computeNormal(depth0, depthIntrinsicsInv0, x, y);
			const float cInput = intensity0(x, y);
			const float gInput = gradmag0(x, y);
			if (p0.x != INVALID && n0.x != INVALID && gInput != INVALID && d0 >= depthMin && d0 <= depthMax) {
				const vec3f pTransInput = transform0to1 * p0;
				const vec3f nTransInput = transform0to1.getRotation() * n0;
				vec2f screenPosf = cameraToDepth(depthIntrinsics1, pTransInput);
				vec2i screenPos = math::round(screenPosf);
#ifdef DEBUG_IMAGES_REPROJ
				pc0trans.m_points.push_back(pTransInput);
				pc0trans.m_colors.push_back(vec4f(cInput, cInput, cInput, 1.0f));
#endif
				if (screenPos.x >= 0 && screenPos.y >= 0 && screenPos.x < (int)depth0.getWidth() && screenPos.y < (int)depth0.getHeight()) {
					const float dTarget = depth1(screenPos.x, screenPos.y);
					const vec3f pTarget = computeCameraSpacePosition(depth1, depthIntrinsicsInv1, screenPos.x, screenPos.y);
					const vec3f nTarget = computeNormal(depth1, depthIntrinsicsInv1, screenPos.x, screenPos.y);
					const float cTarget = intensity1(screenPos.x, screenPos.y);
					const float gTarget = gradmag1(screenPos.x, screenPos.y);

					if (pTarget.x != INVALID && nTarget.x != INVALID && gTarget != INVALID && dTarget >= depthMin && dTarget <= depthMax) {
						float d = (pTransInput - pTarget).length();
						float dNormal = nTransInput | nTarget;
						float g = std::fabs(gInput - gTarget);
#ifdef DEBUG_IMAGES_REPROJ
						pc1.m_points.push_back(pTarget);
						pc1.m_colors.push_back(vec4f(cTarget, cTarget, cTarget, 1.0f));
#endif
						if (dNormal >= normalThresh && d <= distThresh /*&& g <= colorThresh*/) { 
							//const float weight = std::max(0.0f, 0.5f*((1.0f - d / distThresh) + (1.0f - cameraToKinectProjZ(pTransInput.z, depthMin, depthMax)))); // for weighted ICP;
							err.numCorrs++;
							err.depthL1 += std::fabs(pTransInput.z - dTarget); //l1 between projected depth and depth
							err.depthL2 += d;	//residual
							err.intensityL1 += std::fabs(cInput - cTarget);
							err.intensityGradL1 += std::fabs(g);
#ifdef DEBUG_IMAGES_REPROJ
							debug(x/subsampleFactor, y/subsampleFactor) = d;
#endif
						}
					} // projected to valid depth
				} // inside image
			}
		} // x
	} // y

#ifdef DEBUG_IMAGES_REPROJ
	FreeImageWrapper::saveImage("depth0.png", ColorImageR32G32B32(depth0));
	FreeImageWrapper::saveImage("depth1.png", ColorImageR32G32B32(depth1));
	FreeImageWrapper::saveImage("intensity0.png", intensity0);
	FreeImageWrapper::saveImage("intensity1.png", intensity1);
	FreeImageWrapper::saveImage("gradmag0.png", gradmag0);
	FreeImageWrapper::saveImage("gradmag1.png", gradmag1);
	PointCloudIOf::saveToFile("pc0trans.ply", pc0trans);
	PointCloudIOf::saveToFile("pc1.ply", pc1);
	FreeImageWrapper::saveImage("corr.png", ColorImageR32G32B32(debug));
	std::cout << "computed corrs" << std::endl;
	std::cout << "\t#corrs = " << err.numCorrs << std::endl;
	if (err.numCorrs > 0) {
		std::cout << "\tdepth l1 = " << err.depthL1 << std::endl;
		std::cout << "\tdepth l1 = " << err.depthL2 << std::endl;
		std::cout << "\tintensity l1 = " << err.intensityL1 << std::endl;
		std::cout << "\tintensity grad l1 = " << err.intensityGradL1 << std::endl;
	}
	getchar();
#endif
	return err;
}

bool Reprojection::checkValidReprojection(const DepthImage32& depth0, const DepthImage32& depth1, const mat4f& depthIntrinsics0, const mat4f& depthIntrinsics1,
	const mat4f& depthIntrinsicsInv0, const mat4f& depthIntrinsicsInv1, const mat4f& transform0to1, float maxDepth)
{
	const unsigned int subsampleFactor = depth0.getWidth() == 640 ? 16 : 32; //todo fix hack
	const float normalThresh = 0.9f;
	const float distThresh = 0.2f;
	const float depthMin = 0.4f;
	const float depthMax = maxDepth;
	const float INVALID = -std::numeric_limits<float>::infinity();

	unsigned int count = 0;

#pragma omp parallel for
	for (int _y = 0; _y < (int)depth0.getHeight(); _y += subsampleFactor) {
		unsigned int y = (unsigned int)_y;
		for (unsigned int x = 0; x < depth0.getWidth(); x += subsampleFactor) {
			const float d0 = depth0(x, y);
			const vec3f p0 = computeCameraSpacePosition(depth0, depthIntrinsicsInv0, x, y);
			const vec3f n0 = computeNormal(depth0, depthIntrinsicsInv0, x, y);
			if (p0.x != INVALID && n0.x != INVALID && d0 >= depthMin && d0 <= depthMax) {
				const vec3f pTransInput = transform0to1 * p0;
				const vec3f nTransInput = transform0to1.getRotation() * n0;
				vec2f screenPosf = cameraToDepth(depthIntrinsics1, pTransInput);
				vec2i screenPos = math::round(screenPosf);
				if (screenPos.x >= 0 && screenPos.y >= 0 && screenPos.x < (int)depth0.getWidth() && screenPos.y < (int)depth0.getHeight()) {
					const float dTarget = depth1(screenPos.x, screenPos.y);
					const vec3f pTarget = computeCameraSpacePosition(depth1, depthIntrinsicsInv1, screenPos.x, screenPos.y);
					const vec3f nTarget = computeNormal(depth1, depthIntrinsicsInv1, screenPos.x, screenPos.y);

					if (pTarget.x != INVALID && nTarget.x != INVALID && dTarget >= depthMin && dTarget <= depthMax) {
						float d = (pTransInput - pTarget).length();
						float dNormal = nTransInput | nTarget;
						if (dNormal >= normalThresh && d <= distThresh) {
							count++;
							if (count >= 10) return true;
						}
					} // projected to valid depth
				} // inside image
			}
		} // x
	} // y
	return false;
}
