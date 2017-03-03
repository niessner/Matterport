// sensorFile.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iomanip>
#include <valarray>

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

void evaluate(const std::string& matterportPath, const std::string& patchFeatPath, const std::vector<std::string>& scenes, const std::string& logFile = "");

void computeTrainStatistics(const std::string& basePath, const std::string& sceneFileList)
{
	const vec2ui imageBounds = vec2ui(640, 480);
	//const vec2ui imageBounds = vec2ui(640, 512);

	std::vector<std::string> scenes;
	{
		std::ifstream s(sceneFileList);
		MLIB_ASSERT(s.is_open());
		std::string line;
		while (std::getline(s, line))
			scenes.push_back(line);
	}
	std::vector<size_t> numMatchesPerScene; const unsigned int padding = 128 / 2; const unsigned int matchSkip = 4;
	std::vector<std::vector<std::pair<vec2f, vec3f>>> keyposPos, keyposNeg; //<image pos, world pos> per scene
	for (const auto& scene : scenes) {
		const std::string matchesFile = basePath + scene + "/matches.txt";
		const std::string negativFile = basePath + scene + "/negatives.txt";
		if (util::fileExists(matchesFile)) {
			keyposPos.push_back(std::vector<std::pair<vec2f, vec3f>>());
			std::ifstream s(matchesFile); std::string line; size_t lineCount = 0;
			while (std::getline(s, line)) {
				if (lineCount >= 3 && ((lineCount - 3) % matchSkip == 0 || (lineCount - 3) % matchSkip == 1)) {
					const auto parts = util::split(line, '\t');
					const auto partsPixelPos = util::split(parts[3], ' ');
					vec2f pixPos(util::convertTo<float>(partsPixelPos[0]), util::convertTo<float>(partsPixelPos[1]));
					pixPos *= 0.5f;
					const auto partsWorldPos = util::split(parts[5], ' ');
					vec3f worldPos(util::convertTo<float>(partsWorldPos[0]), util::convertTo<float>(partsWorldPos[1]), util::convertTo<float>(partsWorldPos[2]));
					keyposPos.back().push_back(std::make_pair(pixPos, worldPos));
				}
				lineCount++;
			}
		}
		if (util::fileExists(negativFile)) {
			keyposNeg.push_back(std::vector<std::pair<vec2f, vec3f>>());
			std::ifstream s(negativFile); std::string line; size_t lineCount = 0;
			while (std::getline(s, line)) {
				if (lineCount >= 3 && ((lineCount - 3) % matchSkip == 0 || (lineCount - 3) % matchSkip == 1)) {
					const auto parts = util::split(line, '\t');
					const auto partsPixelPos = util::split(parts[3], ' ');
					vec2f pixPos(util::convertTo<float>(partsPixelPos[0]), util::convertTo<float>(partsPixelPos[1]));
					pixPos *= 0.5f;
					const auto partsWorldPos = util::split(parts[5], ' ');
					vec3f worldPos(util::convertTo<float>(partsWorldPos[0]), util::convertTo<float>(partsWorldPos[1]), util::convertTo<float>(partsWorldPos[2]));
					keyposNeg.back().push_back(std::make_pair(pixPos, worldPos));
				}
				lineCount++;
			}
		}
	}
	for (unsigned int p = 0; p < keyposPos.size(); p++) {
		const auto& keysPos = keyposPos[p];
		const auto& keysNeg = keyposNeg[p];
		size_t count = 0;
		for (unsigned int i = 0; i < keysPos.size(); i += 2) {
			const vec2f& pixPosA = keysPos[i].first;
			const vec2f& pixPosB = keysPos[i + 1].first;
			const vec2f& pixNegA = keysNeg[i].first;
			const vec2f& pixNegB = keysNeg[i + 1].first;
			if (pixPosA.x >= padding && pixPosA.y >= padding && pixPosA.x + padding <= imageBounds.x && pixPosA.y + padding <= imageBounds.y &&
				pixPosB.x >= padding && pixPosB.y >= padding && pixPosB.x + padding <= imageBounds.x && pixPosB.y + padding <= imageBounds.y &&
				pixNegA.x >= padding && pixNegA.y >= padding && pixNegA.x + padding <= imageBounds.x && pixNegA.y + padding <= imageBounds.y &&
				pixNegB.x >= padding && pixNegB.y >= padding && pixNegB.x + padding <= imageBounds.x && pixNegB.y + padding <= imageBounds.y) {
				count++; //#valid key matches per scene
				//stats
			}
		}
		numMatchesPerScene.push_back(count);
	}
	std::sort(numMatchesPerScene.begin(), numMatchesPerScene.end());
	size_t total = std::accumulate(numMatchesPerScene.begin(), numMatchesPerScene.end(), 0);
	float avg = (float)total / (float)numMatchesPerScene.size();
	std::cout << numMatchesPerScene << std::endl;
	std::cout << "total: " << total << std::endl;
	std::cout << "avg: " << avg << std::endl;
}

int main(int argc, char* argv[])
{
#if defined(DEBUG) | defined(_DEBUG)
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
	//_CrtSetBreakAlloc(1042);

	std::srand(0);

	try {
		if (false) {
			const std::string basePath = "//tirion/share/datasets/Matterport/Matching/";
			const std::string sceneFileList = "//tirion/share/datasets/Matterport/Matching/scenes_trainval.txt";
			computeTrainStatistics(basePath, sceneFileList);
			std::cout << "DONE" << std::endl;
			getchar();
		}

		const std::string matterportPath = "//tirion/share/datasets/Matterport/Matching/";
		const std::string patchFeatPath = "//tirion/share/adai/code/Matterport/2dmatch/output/output-test-1-100948/";
		const std::string trainTestPrefix = "scenes_";
		const std::vector<std::string> phases = { "test" }; //the current one is train/test for same scene

		for (const std::string& phase : phases) {
			const std::string fileList = matterportPath + trainTestPrefix + phase + ".txt";
			std::vector<std::string> scenes = readScenes(fileList);
			if (scenes.empty()) { std::cout << "warning: no scenes found from " << fileList << std::endl; continue; }
			evaluate(matterportPath, patchFeatPath + phase + "/", scenes, phase + ".csv");
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

//#define USE_POS_NEG
void evaluate(const std::string& matterportPath, const std::string& patchFeatPath, const std::vector<std::string>& scenes, const std::string& logFile /*= ""*/)
{
	//TODO: output stats per scenes?
	//Note: currently ignores scenes

	const unsigned int numDescDistThresholds = 1000;
	Timer t;

	std::vector<float> ancPosDists, ancNegDists, posNegDists;
	if (true) {	//read feature vectors 
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
#ifdef USE_POS_NEG
		computeDistancesL2(posFeats, negFeats, posNegDists);
#endif
		t.stop(); std::cout << "done! (" << t.getElapsedTime() << " s)" << std::endl;
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

	const auto ancPosMaxIt = std::max_element(ancPosDists.begin(), ancPosDists.end());
	const auto ancNegMaxIt = std::max_element(ancNegDists.begin(), ancNegDists.end());
#ifdef USE_POS_NEG
	const auto posNegMaxIt = std::max_element(posNegDists.begin(), posNegDists.end());
	const float maxDist = std::max(*ancPosMaxIt, std::max(*ancNegMaxIt, *posNegMaxIt));
#else
	const float maxDist = std::max(*ancPosMaxIt, *ancNegMaxIt);
#endif
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
		//if (recall[i] > 0.949 && recall[i] < 0.951) {
		if (recall[i] > 0.94 && recall[i] < 0.96) {
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
}

