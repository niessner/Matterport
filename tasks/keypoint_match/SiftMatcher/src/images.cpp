
#include "stdafx.h"

#include "images.h"
#include "globalAppState.h"

#include "mLibFLANN.h"


void Images::loadMatchFile(const std::string& filename, std::vector<KeyPointMatch>& matches,
	std::vector<KeyPoint>& keypoints, std::vector<std::unordered_map<KeyPoint, size_t>>& keymaps) const
{
	std::ifstream s(filename);
	if (!s.is_open()) throw MLIB_EXCEPTION("failed to open file " + filename + " for read");

	const float scale = (float)m_images.front().front().getWidth() / (float)m_keyDetectImageWidth;
	if (scale != (float)m_images.front().front().getHeight() / (float)m_keyDetectImageHeight)
		throw MLIB_EXCEPTION("must have isotropic scale");

	const unsigned int skip = GAS::get().s_keyMatchSkip;

	std::string line;
	size_t lineCount = 0, kpMatchCount = 0;
	while (std::getline(s, line)) {
		if (lineCount >= 3) { //skip header lines
			if ((lineCount - 3) % skip == 0 || (lineCount - 3) % skip == 1) {
				KeyPoint kp; vec2f offset;

				const auto parts = util::split(line, '\t');
				kp.m_sensorIdx = util::convertTo<unsigned int>(parts[1]);
				kp.m_imageIdx = util::convertTo<unsigned int>(parts[2]);
				kp.m_pixelPos = util::convertTo<vec2f>(parts[3]);
				kp.m_depth = util::convertTo<float>(parts[4]);
				kp.m_worldPos = util::convertTo<vec3f>(parts[5]);
				offset = util::convertTo<float>(parts[6]);
				kp.m_size = util::convertTo<float>(parts[7]);
				kp.m_angle = util::convertTo<float>(parts[8]);
				kp.m_octave = util::convertTo<int>(parts[9]);
				kp.m_scale = util::convertTo<float>(parts[10]);
				kp.m_opencvPackOctave = util::convertTo<int>(parts[11]);
				kp.m_response = -std::numeric_limits<float>::infinity(); //not recorded in matches.txt

				kp.m_pixelPos *= scale; 
				kp.m_size *= scale;

				//add to keypoints list
				size_t idx;
				auto& keymap = keymaps[kp.m_sensorIdx];
				const auto it = keymap.find(kp);
				if (it == keymap.end()) {
					idx = keypoints.size();
					keymap[kp] = idx;
					keypoints.push_back(kp);
				}
				else {
					idx = it->second;
				}

				if (kpMatchCount % 2 == 0) { //new match
					matches.push_back(KeyPointMatch());
					matches.back().m_kp0 = idx;
				}
				else { //part of current match
					matches.back().m_kp1 = idx;
					matches.back().m_offset = offset;
				}
				kpMatchCount++;
				std::cout << "\r[ read " << kpMatchCount << " ]";
			}
		}
		lineCount++;
	}
	std::cout << "\r[ read " << kpMatchCount << " of " << lineCount << " lines ]" << std::endl;
}


void Images::loadGTMatches(const std::string& filenamePos, const std::string& filenameNeg)
{
	m_keyPoints.clear();
	m_gtKeyPointMatchesPos.clear();
	m_gtKeyPointMatchesNeg.clear();

	std::vector<std::unordered_map<KeyPoint, size_t>> keymaps(GAS::get().s_maxNumSensors);

	//load positive matches
	loadMatchFile(filenamePos, m_gtKeyPointMatchesPos, m_keyPoints, keymaps);
	//load negative matches
	loadMatchFile(filenameNeg, m_gtKeyPointMatchesNeg, m_keyPoints, keymaps);

	const unsigned int pad = GAS::get().s_keyPadding / 2;
	const unsigned int imageWidth = m_images.front().front().getWidth();
	const unsigned int imageHeight = m_images.front().front().getHeight();

	std::vector<KeyPointMatch> poss, negs;
	for (unsigned int i = 0; i < m_gtKeyPointMatchesPos.size(); i++) {
		const auto& pos = m_gtKeyPointMatchesPos[i];
		const auto& neg = m_gtKeyPointMatchesNeg[i];
		MLIB_ASSERT(pos.m_kp0 == neg.m_kp0);
		const auto& kpAnc = m_keyPoints[pos.m_kp0];
		const auto& kpPos = m_keyPoints[pos.m_kp1];
		const auto& kpNeg = m_keyPoints[pos.m_kp1];

		const bool bInBounds = 
			kpAnc.m_pixelPos.x >= pad &&  kpAnc.m_pixelPos.x + pad <= imageWidth &&
			kpAnc.m_pixelPos.y >= pad &&  kpAnc.m_pixelPos.y + pad <= imageHeight &&
			kpPos.m_pixelPos.x >= pad &&  kpPos.m_pixelPos.x + pad <= imageWidth &&
			kpPos.m_pixelPos.y >= pad &&  kpPos.m_pixelPos.y + pad <= imageHeight &&
			kpNeg.m_pixelPos.x >= pad &&  kpNeg.m_pixelPos.x + pad <= imageWidth &&
			kpNeg.m_pixelPos.y >= pad &&  kpNeg.m_pixelPos.y + pad <= imageHeight;

		if (bInBounds) { 
			poss.push_back(pos);
			negs.push_back(neg);
		}
	}
	m_gtKeyPointMatchesPos = poss;
	m_gtKeyPointMatchesNeg = negs;

	//MatchVisualization::visulizeMatches(m_images, m_gtKeyPointMatchesPos, m_keyPoints, 10);
	//MatchVisualization::visulizeMatches(m_images, m_gtKeyPointMatchesNeg, m_keyPoints, 10);
	//int a = 5;
}


void Images::matchKeyPoints(const std::string& logFile, bool dumpTemporaryBinary, const std::string& dataPath /*= ""*/) const
{
	//const float siftMatchThresh = GAS::get().s_siftMatchThresh;

	Timer t;
	std::cout << "computing feature distances... "; t.start();
	std::vector<float> matchDistsPos, matchDistsNeg;
	KeyPointMatcher::matchKeyPoints(m_images, m_keyPoints, m_gtKeyPointMatchesPos, matchDistsPos, GAS::get().s_featureType);
	KeyPointMatcher::matchKeyPoints(m_images, m_keyPoints, m_gtKeyPointMatchesNeg, matchDistsNeg, GAS::get().s_featureType);
	t.stop(); std::cout << "done! (" << t.getElapsedTime() << " s)" << std::endl;

	if (dumpTemporaryBinary) {
		if (!util::directoryExists(dataPath)) util::makeDirectory(dataPath);
		const std::string file = dataPath + m_name + ".bin";
		BinaryDataStreamFile s(file, true);
		s << matchDistsPos;
		s << matchDistsNeg;
		s.close();
	}
	else {
		computePrecisionRecall(matchDistsPos, matchDistsNeg, logFile);
	}

}

ml::vec4ui Images::evaluatePosMatch(const std::vector<float>& dists, const float thresh)
{
	unsigned int truePos = 0, falsePos = 0, trueNeg = 0, falseNeg = 0;
	for (const float d : dists) {
		if (d < thresh) truePos++;
		else falsePos++;
	}
	return vec4ui(truePos, falsePos, trueNeg, falseNeg);
}

ml::vec4ui Images::evaluateNegMatch(const std::vector<float>& dists, const float thresh)
{
	unsigned int truePos = 0, falsePos = 0, trueNeg = 0, falseNeg = 0;
	for (const float d : dists) {
		if (d < thresh) falseNeg++;
		else trueNeg++;
	}
	return vec4ui(truePos, falsePos, trueNeg, falseNeg);
}

void Images::computeKeyMatchStatistics(const std::string& logFile) const
{
	if (m_gtKeyPointMatchesPos.empty() || m_gtKeyPointMatchesNeg.empty()) {
		std::cout << "no keypoint matches to compute stats for" << std::endl;
		return;
	}
	//const std::string outdir = "stats/";
	//if (!util::directoryExists(outdir)) util::makeDirectory(outdir);

	//PointCloudf pcAnc, pcPos, pcNeg; //point clouds of keypoints involved in matches (there may be duplicates)
	SparseGrid3<unsigned int> grid; // hash all keypoints of matches 
	const float hashVoxSize = 0.1f; //10cm voxels

	const vec4f colorAnc = vec4f(0.0f, 0.0f, 1.0f, 1.0f); //anc -> blue
	const vec4f colorPos = vec4f(0.0f, 1.0f, 0.0f, 1.0f); //pos -> green
	const vec4f colorNeg = vec4f(1.0f, 0.0f, 0.0f, 1.0f); //neg -> red

	//compute bbox
	bbox3f bbox;
	for (unsigned int i = 0; i < m_gtKeyPointMatchesPos.size(); i++) {
		const KeyPointMatch& pos = m_gtKeyPointMatchesPos[i];
		const KeyPointMatch& neg = m_gtKeyPointMatchesNeg[i];
		const vec3f& worldPosAnc = m_keyPoints[pos.m_kp0].m_worldPos;
		const vec3f& worldPosPos = m_keyPoints[pos.m_kp1].m_worldPos;
		const vec3f& worldPosNeg = m_keyPoints[neg.m_kp1].m_worldPos;
		bbox.include(worldPosAnc);
		bbox.include(worldPosPos);
		bbox.include(worldPosNeg);
		//bbox.include(m_keyPoints[m_gtKeyPointMatchesPos[i].m_kp0].m_worldPos);
		//bbox.include(m_keyPoints[m_gtKeyPointMatchesPos[i].m_kp1].m_worldPos);
		//bbox.include(m_keyPoints[m_gtKeyPointMatchesNeg[i].m_kp1].m_worldPos);
	}
	const vec3i gridDim = math::round(bbox.getExtent() / hashVoxSize);
	mat4f worldToVoxel = bbox.worldToCubeTransform(); //[-1,1]
	worldToVoxel = mat4f::scale(0.5f*vec3f(gridDim - 1)) * mat4f::translation(1.0f) * worldToVoxel; //[0,1]

	std::unordered_map<unsigned int, unsigned int> numMatchesPerKey;
	for (unsigned int i = 0; i < m_gtKeyPointMatchesPos.size(); i++) {
		const KeyPointMatch& pos = m_gtKeyPointMatchesPos[i];
		const KeyPointMatch& neg = m_gtKeyPointMatchesNeg[i];

		const vec3f& worldPosAnc = m_keyPoints[pos.m_kp0].m_worldPos;
		const vec3f& worldPosPos = m_keyPoints[pos.m_kp1].m_worldPos;
		const vec3f& worldPosNeg = m_keyPoints[neg.m_kp1].m_worldPos;

		//pcAnc.m_points.push_back(worldPosAnc);	pcAnc.m_colors.push_back(colorAnc);
		//pcPos.m_points.push_back(worldPosPos);	pcPos.m_colors.push_back(colorPos);
		//pcNeg.m_points.push_back(worldPosNeg);	pcNeg.m_colors.push_back(colorNeg);

		const vec3i voxPosAnc = math::round(worldToVoxel * worldPosAnc);
		const vec3i voxPosPos = math::round(worldToVoxel * worldPosPos);
		const vec3i voxPosNeg = math::round(worldToVoxel * worldPosNeg);
		if (grid.exists(voxPosAnc)) grid[voxPosAnc]++; else grid[voxPosAnc] = 1;
		if (grid.exists(voxPosPos)) grid[voxPosPos]++; else grid[voxPosPos] = 1;
		if (grid.exists(voxPosNeg)) grid[voxPosNeg]++; else grid[voxPosNeg] = 1;

		auto it = numMatchesPerKey.find(pos.m_kp0);
		if (it == numMatchesPerKey.end()) numMatchesPerKey[pos.m_kp0] = 1;
		else it->second++;
	}

	unsigned int numEntries = 0, numVoxels = 0;
	for (auto it = grid.begin(); it != grid.end(); ++it) {
		numVoxels++;
		numEntries += it->second;
	}
	std::cout << "#points/voxel (" << hashVoxSize << " m): " << (float)numEntries / (float)numVoxels << " = (" << numEntries << " / " << numVoxels << ")" << std::endl;
	std::cout << "\tgrid size = " << gridDim << " (" << gridDim.x*gridDim.y*gridDim.z << ")" << std::endl;

	std::vector<vec2ui> numMatchesPerKeyVec; numMatchesPerKeyVec.reserve(numMatchesPerKey.size());
	float avgNumMatchesPerKey = 0.0f;
	for (const auto m : numMatchesPerKey) {
		numMatchesPerKeyVec.push_back(vec2ui(m.first, m.second));
		avgNumMatchesPerKey += m.second;
	}
	avgNumMatchesPerKey /= (float)numMatchesPerKey.size();
	std::sort(numMatchesPerKeyVec.begin(), numMatchesPerKeyVec.end(), [](const vec2ui &left, const vec2ui &right) {
		return left.y > right.y;
	});
	std::cout << "max #matches per key = " << numMatchesPerKeyVec.front() << std::endl;

	//PointCloudIOf::saveToFile(outdir + "anc.ply", pcAnc);
	//PointCloudIOf::saveToFile(outdir + "pos.ply", pcPos);
	//PointCloudIOf::saveToFile(outdir + "neg.ply", pcNeg);

	const bool bPrintHeader = !util::fileExists(logFile);
	std::ofstream ofs(logFile, std::ios::app); const std::string splitter = ",";
	if (!ofs.is_open()) throw MLIB_EXCEPTION("failed to open " + logFile + " for write");
	if (bPrintHeader) ofs << "scene_name" << splitter << "grid dim (vox size: " << hashVoxSize << "m)" << splitter << "#points/voxel" << splitter << "#points" << splitter << "#voxels" << splitter << "max #matches per key" << splitter << "avg #matches per key" << std::endl;
	ofs << m_name << splitter << gridDim << splitter << (float)numEntries / (float)numVoxels << splitter << numEntries << splitter << numVoxels << splitter << numMatchesPerKeyVec.front().y << splitter << avgNumMatchesPerKey << std::endl;
	ofs.close();
}

void Images::computePrecisionRecall(const std::vector<float> &matchDistsPos, const std::vector<float> &matchDistsNeg, const std::string& logFile)
{
	const unsigned int numThresholds = GAS::get().s_numThresholds;

	Timer t;
	//create thresholds
	const auto maxDistPosIt = std::max_element(matchDistsPos.begin(), matchDistsPos.end());
	const auto maxDistNegIt = std::max_element(matchDistsNeg.begin(), matchDistsNeg.end());
	const float maxDist = std::max(*maxDistPosIt, *maxDistNegIt);
	std::cout << "==> max dist = " << maxDist << std::endl;
	std::vector<float> thresholds(numThresholds);
	const float inc = 1.0f / (float)(numThresholds - 1); //thresholds range [0,maxDist]
	for (unsigned int i = 0; i < numThresholds; i++) thresholds[i] = (float)i * inc * maxDist;

	std::cout << "computing precision/recall... "; t.start();
	std::vector<vec4ui> results(thresholds.size()); // vec4<true positive, false positive, true negative, false negative>
	std::vector<float> accuracy(thresholds.size()), precision(thresholds.size()), recall(thresholds.size()), rateFP(thresholds.size());
#pragma omp parallel for
	for (int i = 0; i < (int)thresholds.size(); i++) {
		vec4ui resAncPos = evaluatePosMatch(matchDistsPos, thresholds[i]);
		vec4ui resAncNeg = evaluateNegMatch(matchDistsNeg, thresholds[i]);

		results[i] = resAncPos + resAncNeg;
		accuracy[i] = (float)(results[i][0] + results[i][2]) / (float)(results[i][0] + results[i][1] + results[i][2] + results[i][3]);
		precision[i] = (float)results[i][0] / (float)(results[i][0] + results[i][1]);
		recall[i] = (float)results[i][0] / (float)(results[i][0] + results[i][3]);
		rateFP[i] = (float)results[i][1] / (float)(results[i][1] + results[i][2]);
	}
	float errAt95Recall = 0.0f; unsigned int norm = 0;
	for (unsigned int i = 0; i < thresholds.size(); i++) {
		if (recall[i] > 0.949 && recall[i] < 0.951) {
		//if (recall[i] > 0.94 && recall[i] < 0.96) {
			errAt95Recall += rateFP[i];
			norm++;
		}
	}
	errAt95Recall /= (float)norm;
	t.stop(); std::cout << "done! (" << t.getElapsedTime() << " s)" << std::endl;
	std::cout << "False-positive rate (error) at 95% recall: " << errAt95Recall << std::endl;

	//log dists to file
	std::cout << "logging to file: " << logFile << "... ";
	std::ofstream s(logFile, std::ios::app); const std::string splitter = ",";
	if (!s.is_open()) throw MLIB_EXCEPTION("failed to open " + logFile + " for write");
	s << "threshold" << splitter << "#true positive" << splitter << "#false positive" << splitter << "#true negative" << splitter << "#false negative" << splitter << "accuracy" << splitter << "precision" << splitter << "recall" << std::endl;
	for (unsigned int i = 0; i < thresholds.size(); i++) {
		s << thresholds[i] << splitter << results[i][0] << splitter << results[i][1] << splitter << results[i][2] << splitter << results[i][3] << splitter << accuracy[i] << splitter << precision[i] << splitter << recall[i] << std::endl;
	}
	std::cout << "done!" << std::endl;
}

void Images::computePrecisionRecall(const std::string& logFile, const std::string& dataPath /*= ""*/)
{
	if (!util::directoryExists(dataPath)) {
		std::cout << "data dir (" << dataPath << ") does not exist!" << std::endl;
		return;
	}
	std::vector<float> distsPos, distsNeg;

	Directory dir(dataPath);
	const auto files = dir.getFiles();
	for (const auto& f : files) {
		const std::string file = dataPath + f;
		if (util::fileExists(file)) {
			std::vector<float> pos, neg;
			BinaryDataStreamFile s(file, false);
			s >> pos;
			s >> neg;
			s.close();
			distsPos.insert(distsPos.end(), pos.begin(), pos.end());
			distsNeg.insert(distsNeg.end(), neg.begin(), neg.end());
		}
	}
	computePrecisionRecall(distsPos, distsNeg, logFile);
}
