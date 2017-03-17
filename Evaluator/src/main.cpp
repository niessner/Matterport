// sensorFile.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iomanip>
#include <valarray>
enum DISTANCE_TYPE {
	L2_DIST
};

void generateTrainTestSplit(const std::string matterportPath);

std::vector<std::string> readScenes(const std::string& filename)
{
	if (!util::fileExists(filename)) throw MLIB_EXCEPTION("file " + filename + " does not exist!");

	std::ifstream s(filename);
	if (!s.is_open()) throw MLIB_EXCEPTION("failed to open " + filename + " for read");
	std::vector<std::string> scenes;
	std::string line;
	while (std::getline(s, line)) {
		scenes.push_back(line);
	}
	return scenes;
}

float evaluate(const std::string& matterportPath, const std::string& patchFeatPath, const std::vector<std::string>& scenes,
	DISTANCE_TYPE distType, const std::vector<float>& camAnglesPos, const std::vector<float>& camAnglesNeg,
	const float angleThresh = std::numeric_limits<float>::infinity(), const std::string& logFile = "");
void loadCameras(const std::vector<std::string>& scenes, std::unordered_map< std::string, std::vector<std::vector<mat4f>> >& trajectories);
void computeCameraAngles(const std::string& matchesFile, const std::string& negsFile,
	const std::vector<std::vector<mat4f>>& trajectories,
	std::vector<float>& cameraAnglesPos, std::vector<float>& cameraAnglesNeg);

int main(int argc, char* argv[])
{
#if defined(DEBUG) | defined(_DEBUG)
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
	//_CrtSetBreakAlloc(1042);

	std::srand(0);

	if (false) {
		const bool bCountOnly = true;
		const std::string srcPath = "//tirion/share/datasets/Matterport/Matching1/";
		const std::string dstPath = "//moonraker/adai/data/Matterport/Matching1/";
		//const std::string dstPath = "//doriath/share/datasets/Matterport/Matching1/";
		const std::vector<std::string> files = { "matches1.txt", "negatives1.txt" };
		Directory dirSrc(srcPath);
		const auto& scenes = dirSrc.getDirectories();
		unsigned int count = 0;
		std::vector<std::string> missingImagesScenes;
		for (const auto& scene : scenes) {
			if (bCountOnly) {
				bool done = true;
				for (const auto& f : files) {
					if (!util::fileExists(srcPath + scene + "/" + f)) {
						done = false;
						std::cout << "missing: " << scene << std::endl;
						break;
					}
				}
				if (done) count++;
			}
			else {
				if (!util::directoryExists(dstPath + scene))
					util::makeDirectory(dstPath + scene);
				for (const auto& f : files) {
					if (util::fileExists(srcPath + scene + "/" + f)) {
						util::copyFile(srcPath + scene + "/" + f, dstPath + scene + "/" + f);
						count++;
						std::cout << "\r" << count << "\t" << scene;
					}
				}
				if (!util::directoryExists(dstPath + scene + "/images") || (util::getFileSize(dstPath + scene + "/images/color-00-000000.jpg") != util::getFileSize(srcPath + scene + "/images/color-00-000000.jpg"))) {
					missingImagesScenes.push_back(scene);
				}
			}
		}
		if (bCountOnly)
			std::cout << "count = " << count << std::endl;
		std::cout << std::endl;
		std::cout << "#missing images scenes = " << missingImagesScenes.size() << std::endl;
		if (!missingImagesScenes.empty()) std::cout << missingImagesScenes << std::endl;
		std::cout << "done" << std::endl;
		getchar();
	}

	try {
		const bool useCameras = true;
		const std::string matterportPath = "//tirion/share/datasets/Matterport/MatchingSun3d/";
		//const std::string patchFeatPath = "//tirion/share/adai/code/Matterport/2dmatch/output-sun3d/all_1-2000/"; const DISTANCE_TYPE distType = DISTANCE_TYPE::L2_DIST;
		const std::string patchFeatPath = "//moonraker/adai/code/Matterport/2dmatch/output-mp/s3/"; const DISTANCE_TYPE distType = DISTANCE_TYPE::L2_DIST;
		//const std::string patchFeatPath = "//doriath/share/adai/code/Matterport/2dmatch/output-sun3d/mps3-ms800_1-6000/";     const DISTANCE_TYPE distType = DISTANCE_TYPE::L2_DIST;
		//const std::string patchFeatPath = "//tirion/share/adai/code/Matterport/matchnet/output/sun3d/tmp/"; const DISTANCE_TYPE distType = DISTANCE_TYPE::L2_DIST;
		const std::string trainTestPrefix = "scenes_";
		const std::vector<std::string> phases = { "test" }; //the current one is train/test for same scene

		for (const std::string& phase : phases) {
			const std::string fileList = matterportPath + trainTestPrefix + phase + ".txt";
			std::vector<std::string> scenes = readScenes(fileList);
			if (scenes.empty()) { std::cout << "warning: no scenes found from " << fileList << std::endl; continue; }
			std::unordered_map< std::string, std::vector<std::vector<mat4f>> > trajectories;
			std::vector<float> camAnglesPos, camAnglesNeg;
			if (useCameras) {
				const std::string camAngleCacheFile = "cameraAngles.cache";
				if (util::fileExists(camAngleCacheFile)) {
					BinaryDataStreamFile s(camAngleCacheFile, false);
					s >> camAnglesPos; s >> camAnglesNeg; s.close();
				}
				else {
					loadCameras(scenes, trajectories);
					for (const auto& a : trajectories) {
						computeCameraAngles(matterportPath + a.first + "/matches1.txt", matterportPath + a.first + "/negatives1.txt",
							a.second, camAnglesPos, camAnglesNeg);
					}
					BinaryDataStreamFile s(camAngleCacheFile, true);
					s << camAnglesPos; s << camAnglesNeg; s.close();
				}
				MLIB_ASSERT(camAnglesPos.size() == camAnglesNeg.size());
				const unsigned int numAngleThresh = 100;
				const float maxAngle = std::max(*std::max_element(camAnglesPos.begin(), camAnglesPos.end()), *std::max_element(camAnglesNeg.begin(), camAnglesNeg.end()));
				std::cout << "max angle = " << maxAngle << std::endl;
				std::ofstream ofs("angles.csv"); ofs << "angle thresh (degrees),fp" << std::endl;
				for (unsigned int i = 1; i <= numAngleThresh; i++) {
					const float angleThresh = maxAngle*(float)i/(float)numAngleThresh;
					const float fp = evaluate(matterportPath, patchFeatPath + phase + "/", scenes, distType, camAnglesPos, camAnglesNeg, angleThresh);
					if (fp != 0) ofs << math::radiansToDegrees(angleThresh) << "," << fp << std::endl;
				}
				ofs.close();
				std::cout << "DONE DONE DONE" << std::endl;
			}
			else {
				evaluate(matterportPath, patchFeatPath + phase + "/", scenes, distType, camAnglesPos, camAnglesNeg,
					std::numeric_limits<float>::infinity(), phase + ".csv");
			}
		}


		//const std::string matterportPath = "W:/data/matterport/v2/";
		//const bool bGenTrainTestSplit = false;
		//if (bGenTrainTestSplit) { //generate train/test splits
		//	generateTrainTestSplit(matterportPath);
		//}
	}
	catch (const std::exception& e)
	{
		MessageBoxA(NULL, e.what(), "Exception caught", MB_ICONERROR);
		exit(EXIT_FAILURE);
	}
	catch (...)
	{
		MessageBoxA(NULL, "UNKNOWN EXCEPTION", "Exception caught", MB_ICONERROR);
		exit(EXIT_FAILURE);
	}

	std::cout << "<press key to continue>" << std::endl;
	getchar();
	return 0;
}

void generateTrainTestSplit(const std::string matterportPath)
{
	const float percentTrain = 0.7f; const float percentVal = 0.1f;// const float percentTest = 0.2f;
	const std::string outTrainFile = "scenes_train.txt";
	const std::string outValFile = "scenes_val.txt";
	const std::string outTestFile = "scenes_test.txt";
	Directory dir(matterportPath);
	std::vector<std::string> scans = dir.getDirectories();
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(scans.begin(), scans.end(), g);
	const unsigned int numTrainScenes = (unsigned int)std::round(percentTrain * (float)scans.size());
	const unsigned int numValScenes = (unsigned int)std::round(percentVal * (float)scans.size());
	const unsigned int numTestScenes = (unsigned int)scans.size() - numTrainScenes - numValScenes;
	if (numTrainScenes > 0) {
		std::ofstream s(outTrainFile);
		for (unsigned int i = 0; i < numTrainScenes; i++)
			s << scans[i] << std::endl;
	}
	if (numValScenes > 0) {
		std::ofstream s(outValFile);
		for (unsigned int i = numTrainScenes; i < numTrainScenes + numValScenes; i++)
			s << scans[i] << std::endl;
	}
	if (numTestScenes > 0) {
		std::ofstream s(outTestFile);
		for (unsigned int i = numTrainScenes + numValScenes; i < numTrainScenes + numValScenes + numTestScenes; i++)
			s << scans[i] << std::endl;
	}
	std::cout << "generated " << numTrainScenes << " train, " << numValScenes << " val scenes, " << numTestScenes << " test scans" << std::endl;
}

struct KeyMatch {
	vec2ui sensImIdx0;
	vec2ui sensImIdx1;
	vec2f pixPos0;
	vec2f pixPos1;
	vec3f worldPos0;
	vec3f worldPos1;
};
void readKeyMatches(const std::string& matchesFile, std::vector<KeyMatch>& res)
{
	const float scale = 1.0f;
	const unsigned int skip = 1;

	std::ifstream s(matchesFile);
	if (!s.is_open()) throw MLIB_EXCEPTION("failed to open file " + matchesFile + " for read");
	std::string line;
	size_t lineCount = 0, kpMatchCount = 0;
	while (std::getline(s, line)) {
		if (lineCount >= 3) { //skip header lines
			if ((lineCount - 3) % skip == 0 || (lineCount - 3) % skip == 1) {
				const auto parts = util::split(line, '\t');
				const unsigned int sensorIdx = util::convertTo<unsigned int>(parts[1]);
				const unsigned int imageIdx = util::convertTo<unsigned int>(parts[2]);
				vec2f pixelPos = util::convertTo<vec2f>(parts[3]);
				const vec3f worldPos = util::convertTo<vec3f>(parts[5]);
				pixelPos *= scale;

				if (kpMatchCount % 2 == 0) { //new match
					res.push_back(KeyMatch());
					res.back().sensImIdx0 = vec2ui(sensorIdx, imageIdx);
					res.back().pixPos0 = pixelPos;
					res.back().worldPos0 = worldPos;
				}
				else { //part of current match
					res.back().sensImIdx1 = vec2ui(sensorIdx, imageIdx);
					res.back().pixPos1 = pixelPos;
					res.back().worldPos1 = worldPos;
				}
				kpMatchCount++;
				if (kpMatchCount % 1000 == 0) std::cout << "\r[ read " << kpMatchCount << " ]";
			}
		}
		lineCount++;
	}
	std::cout << "\r[ read " << kpMatchCount << " of " << lineCount << " lines ]" << std::endl;
}

void computeCameraAngles(const std::string& matchesFile, const std::string& negsFile,
	const std::vector<std::vector<mat4f>>& trajectories,
	std::vector<float>& cameraAnglesPos, std::vector<float>& cameraAnglesNeg)
{

	std::vector<KeyMatch> keyMatches, keyNegs;
	readKeyMatches(matchesFile, keyMatches);
	readKeyMatches(negsFile, keyNegs);
	std::cout << "#read poss = " << keyMatches.size() << std::endl;
	std::cout << "#read negs = " << keyNegs.size() << std::endl;

	const unsigned int pad = 64 / 2;
	const unsigned int imageWidth = 640;
	const unsigned int imageHeight = 480;

	std::vector<KeyMatch> poss, negs;
	for (unsigned int i = 0; i < keyMatches.size(); i++) {
		const auto& pos = keyMatches[i];
		const auto& neg = keyNegs[i];
		MLIB_ASSERT(pos.sensImIdx0 == neg.sensImIdx0);

		const bool bInBounds =
			pos.pixPos0.x >= pad &&  pos.pixPos0.x + pad <= imageWidth &&
			pos.pixPos0.y >= pad &&  pos.pixPos0.y + pad <= imageHeight &&
			pos.pixPos1.x >= pad &&  pos.pixPos1.x + pad <= imageWidth &&
			pos.pixPos1.y >= pad &&  pos.pixPos1.y + pad <= imageHeight &&
			neg.pixPos1.x >= pad &&  neg.pixPos1.x + pad <= imageWidth &&
			neg.pixPos1.y >= pad &&  neg.pixPos1.y + pad <= imageHeight;

		if (bInBounds) {
			poss.push_back(pos);
			negs.push_back(neg);
		}
	}
	std::cout << "#poss = " << poss.size() << std::endl;
	std::cout << "#negs = " << negs.size() << std::endl;
	//compute angles
	for (unsigned int i = 0; i < poss.size(); i++) {
		const float dotPos = (trajectories[poss[i].sensImIdx0.x][poss[i].sensImIdx0.y].getTranslation() - poss[i].worldPos0).getNormalized() |
			(trajectories[poss[i].sensImIdx1.x][poss[i].sensImIdx1.y].getTranslation() - poss[i].worldPos1).getNormalized();
		const float dotNeg = (trajectories[negs[i].sensImIdx0.x][negs[i].sensImIdx0.y].getTranslation() - negs[i].worldPos0).getNormalized() |
			(trajectories[negs[i].sensImIdx1.x][negs[i].sensImIdx1.y].getTranslation() - negs[i].worldPos1).getNormalized();
		cameraAnglesPos.push_back(std::acos(math::clamp(dotPos, -1.0f, 1.0f)));
		cameraAnglesNeg.push_back(std::acos(math::clamp(dotNeg, -1.0f, 1.0f)));
		if (cameraAnglesPos.back() != cameraAnglesPos.back() || cameraAnglesNeg.back() != cameraAnglesNeg.back()) {
			std::cout << "ERROR invalid camera angle: (" << dotPos << ", " << cameraAnglesPos.back() << ") (" << dotNeg << ", " << cameraAnglesNeg.back() << ")" << std::endl;
			getchar();
		}
	}
	std::cout << "#cpos = " << cameraAnglesPos.size() << std::endl;
	std::cout << "#cneg = " << cameraAnglesNeg.size() << std::endl;
}

void loadCameras(const std::vector<std::string>& scenes, std::unordered_map< std::string, std::vector<std::vector<mat4f>> >& trajectories)
{
	trajectories.clear();
	const std::string cacheFile = "cameras.cache";
	if (util::fileExists(cacheFile)) {
		std::cout << "loading cameras from cache... ";
		BinaryDataStreamFile s(cacheFile, false);
		s >> trajectories; s.close();
		std::cout << "done!" << std::endl;
	}
	else {
		std::cout << "reading cameras from sens files..." << std::endl;
		//const std::string datapath = "W:/data/matterport/v1_converted/";
		const std::string datapath = "W:/data/sun3d/";
		unsigned int idx = 0;
		for (const std::string& scene : scenes) {
			//std::vector<std::vector<mat4f>> trajectory(3);
			//for (unsigned int i = 0; i < 3; i++) {
			//	const std::string sensFile = datapath + scene + "/" + scene + "_" + std::to_string(i) + ".sens";
			//	if (!util::fileExists(sensFile)) throw MLIB_EXCEPTION("sens file (" + sensFile + ") does not exist!");
			//	SensorData sd(sensFile);
			//	trajectory[i].resize(sd.m_frames.size());
			//	for (unsigned int f = 0; f < sd.m_frames.size(); f++)
			//		trajectory[i][f] = sd.m_frames[f].getCameraToWorld();
			//}
			std::vector<std::vector<mat4f>> trajectory(1);
			{
				const std::string sensFile = datapath + scene + "/" + scene + "_ground_truth.sens";
				if (!util::fileExists(sensFile)) throw MLIB_EXCEPTION("sens file (" + sensFile + ") does not exist!");
				SensorData sd(sensFile);
				trajectory[0].resize(sd.m_frames.size());
				for (unsigned int f = 0; f < sd.m_frames.size(); f++)
					trajectory[0][f] = sd.m_frames[f].getCameraToWorld();
			}
			trajectories[scene] = trajectory;
			std::cout << "\r[ " << ++idx << " | " << scenes.size() << " ]";
		}
		BinaryDataStreamFile s(cacheFile, true);
		s << trajectories; s.close();
		std::cout << std::endl << "done!" << std::endl;

	}
}

void readFeatureVectors(const std::string& filename, std::vector<std::valarray<float>>& features)
{
	if (!util::fileExists(filename)) throw MLIB_EXCEPTION("file " + filename + " does not exist!");
	features.clear();

	std::ifstream s(filename);
	if (!s.is_open()) throw MLIB_EXCEPTION("failed to open " + filename + " for read");
	unsigned int featSize = 0;
	std::string line;
	while (std::getline(s, line)) {
		const auto parts = util::split(line, ',');
		if (featSize == 0) featSize = (unsigned int)parts.size();
		else if (parts.size() != featSize) {
			std::cerr << "warning: corrupt feature file at feature " << features.size() << " (" << parts.size() << " vs " << featSize << "), truncating read" << std::endl;
			break;
		}
		features.push_back(std::valarray<float>(featSize));
		auto& feat = features.back();
		for (unsigned int i = 0; i < featSize; i++) {
			feat[i] = util::convertTo<float>(parts[i]);
		}
		if (features.size() % 1000 == 0) std::cout << "\r[ feat " << features.size() << " ]";
	}
	std::cout << "\r[ feat " << features.size() << " ]";
	std::cout << std::endl;
}

float computeDistanceL2Sq(const std::valarray<float>& v0, const std::valarray<float>& v1)
{
	MLIB_ASSERT(v0.size() == v1.size());
	auto diff = v0 - v1;
	return (diff * diff).sum();
}
float computeDistanceL2(const std::valarray<float>& v0, const std::valarray<float>& v1)
{
	float d = computeDistanceL2Sq(v0, v1);
	return std::sqrt(d);
}

// compute feats0.size() l2 dists between corresponding feats of feats0 and feats1
void computeDistancesL2(const std::vector<std::valarray<float>>& feats0, const std::vector<std::valarray<float>>& feats1, std::vector<float>& dists)
{
	MLIB_ASSERT(feats0.size() == feats1.size());
	dists.resize(feats0.size());
#pragma omp parallel for
	for (int i = 0; i < (int)dists.size(); i++)
		dists[i] = computeDistanceL2(feats0[i], feats1[i]);
}

//// compute feats0.size() dists between corresponding feats of feats0 and feats1
//void computeDistances(DISTANCE_TYPE distType, const std::vector<std::valarray<float>>& feats0, const std::vector<std::valarray<float>>& feats1, std::vector<float>& dists)
//{
//	switch (distType)
//	{
//		case DISTANCE_TYPE::L2_DIST:
//			computeDistancesL2(feats0, feats1, dists);
//			break;
//		default:
//			//should never get here
//			break;
//	}
//}


vec4ui evaluatePosMatch(const std::vector<float>& dists, const float thresh)
{
	unsigned int truePos = 0, falsePos = 0, trueNeg = 0, falseNeg = 0;
	for (const float d : dists) {
		if (d < thresh) truePos++;
		else falsePos++;
	}
	return vec4ui(truePos, falsePos, trueNeg, falseNeg);
}
vec4ui evaluateNegMatch(const std::vector<float>& dists, const float thresh)
{
	unsigned int truePos = 0, falsePos = 0, trueNeg = 0, falseNeg = 0;
	for (const float d : dists) {
		if (d < thresh) falseNeg++;
		else trueNeg++;
	}
	return vec4ui(truePos, falsePos, trueNeg, falseNeg);
}


float evaluate(const std::string& matterportPath, const std::string& patchFeatPath, const std::vector<std::string>& scenes,
	DISTANCE_TYPE distType, const std::vector<float>& camAnglesPos, const std::vector<float>& camAnglesNeg, 
	const float angleThresh /*= std::numeric_limits<float>::infinity()*/, const std::string& logFile /*= ""*/)
{
	//TODO: output stats per scenes?
	//Note: currently ignores scenes

	const unsigned int numDescDistThresholds = 1000;
	Timer t;

	std::vector<float> ancPosDists, ancNegDists; //hack...
	if (true) {	//read feature vectors 
		const std::string cachefile = "features.cache";
		if (util::fileExists(cachefile)) {
			BinaryDataStreamFile s(cachefile, false);
			s >> ancPosDists;
			s >> ancNegDists;
		}
		else {
			std::cout << "reading feature vectors..." << std::endl; t.start();
			const std::string ancFeatFile = patchFeatPath + "anc_feat.txt";
			const std::string posFeatFile = patchFeatPath + "pos_feat.txt";
			const std::string negFeatFile = patchFeatPath + "neg_feat.txt";
			std::vector<std::valarray<float>> ancFeats, posFeats, negFeats;
			readFeatureVectors(ancFeatFile, ancFeats);
			readFeatureVectors(posFeatFile, posFeats);
			readFeatureVectors(negFeatFile, negFeats);
			size_t numFeats = std::min(ancFeats.size(), std::min(posFeats.size(), negFeats.size()));
			t.stop(); std::cout << "[ " << t.getElapsedTime() << " s ] read " << numFeats << " features (" << ancFeats.size() << ", " << posFeats.size() << ", " << negFeats.size() << ")" << std::endl;

			//anc-pos -> match, anc-neg -> not match, pos-neg -> not match
			std::cout << "computing feature distances... "; t.start();
			computeDistancesL2(ancFeats, posFeats, ancPosDists);
			computeDistancesL2(ancFeats, negFeats, ancNegDists);
			//compute loss
			float loss = 0.0f;
			for (unsigned int i = 0; i < ancPosDists.size(); i++) {
				loss += -std::log(std::exp(-ancPosDists[i]) / (std::exp(-ancPosDists[i]) + std::exp(-ancNegDists[i])));
			}
			std::cout << "loss = " << (loss / (float)ancPosDists.size()) << " (" << loss << " / " << ancPosDists.size() << ")" << std::endl;
			t.stop(); std::cout << "done! (" << t.getElapsedTime() << " s)" << std::endl;
			BinaryDataStreamFile s(cachefile, true);
			s << ancPosDists; s << ancNegDists;
			s.close();
			std::cout << "CACHED" << std::endl;
			getchar();
		}
	}
	else { //generate random distances
		const unsigned int numDists = 170659;
		const float maxDist = 1.41f;
		ancPosDists.resize(numDists);
		ancNegDists.resize(numDists);
		for (unsigned int i = 0; i < numDists; i++) {
			ancPosDists[i] = math::randomUniform(0.0f, maxDist);
			ancNegDists[i] = math::randomUniform(0.0f, maxDist);
		}
	}
	if (angleThresh != std::numeric_limits<float>::infinity()) {
		for (int i = (int)ancPosDists.size()-1; i >= 0; i--) {
			if (camAnglesPos[i] > angleThresh) ancPosDists.erase(ancPosDists.begin() + i);
			if (camAnglesNeg[i] > angleThresh) ancNegDists.erase(ancNegDists.begin() + i);
		}
		std::cout << "[ " << ancPosDists.size() << ", " << ancNegDists.size() << " ] pos/neg dists (thresh " << angleThresh << ")" << std::endl;
		if (ancPosDists.empty() || ancNegDists.empty()) {
			return 0.0f;
		}
	}

	const auto ancPosMaxIt = std::max_element(ancPosDists.begin(), ancPosDists.end());
	const auto ancNegMaxIt = std::max_element(ancNegDists.begin(), ancNegDists.end());
	const float maxDist = std::max(*ancPosMaxIt, *ancNegMaxIt);
	std::cout << "==> max dist = " << maxDist << std::endl;

	std::vector<float> thresholds(numDescDistThresholds);
	const float inc = 1.0f / (float)(numDescDistThresholds - 1); //thresholds range [0,maxDist]
	for (unsigned int i = 0; i < numDescDistThresholds; i++) thresholds[i] = (float)i * inc * maxDist;

	std::cout << "computing precision/recall... "; t.start();
	std::vector<vec4ui> results(thresholds.size()); // vec4<true positive, false positive, true negative, false negative>
	std::vector<float> accuracy(thresholds.size()), precision(thresholds.size()), recall(thresholds.size()), rateFP(thresholds.size());
#pragma omp parallel for
	for (int i = 0; i < (int)thresholds.size(); i++) {
		vec4ui resAncPos = evaluatePosMatch(ancPosDists, thresholds[i]);
		vec4ui resAncNeg = evaluateNegMatch(ancNegDists, thresholds[i]);
#ifdef USE_POS_NEG
		vec4ui resPosNeg = evaluateNegMatch(posNegDists, thresholds[i]);
		results[i] = resAncPos + resAncNeg + resPosNeg;
#else
		results[i] = resAncPos + resAncNeg;
#endif
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
	t.stop(); std::cout << "done! (" << t.getElapsedTime() << " s)";

	std::cout << std::endl;
	std::cout << "False-positive rate (error) at 95% recall: " << errAt95Recall << std::endl;
	std::cout << std::endl;

	if (!logFile.empty()) {
		std::cout << "logging to file: " << logFile << "... ";
		std::ofstream s(logFile); const std::string splitter = ",";
		if (!s.is_open()) throw MLIB_EXCEPTION("failed to open " + logFile + " for write");
		s << "threshold" << splitter << "#true positive" << splitter << "#false positive" << splitter << "#true negative" << splitter << "#false negative" << splitter << "accuracy" << splitter << "precision" << splitter << "recall" << std::endl;
		for (unsigned int i = 0; i < thresholds.size(); i++) {
			s << thresholds[i] << splitter << results[i][0] << splitter << results[i][1] << splitter << results[i][2] << splitter << results[i][3] << splitter << accuracy[i] << splitter << precision[i] << splitter << recall[i] << std::endl;
		}
		std::cout << "done!" << std::endl;
	}
	return errAt95Recall;
}

