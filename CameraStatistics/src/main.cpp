// sensorFile.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iomanip>
#include <valarray>

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

void loadCameras(const std::string& dataPath, const std::vector<std::string>& scenes, std::unordered_map< std::string, std::vector<std::vector<mat4f>> >& trajectories);
void computeCameraBaselines(const std::string& dataPath, const std::string& matchesFile, const std::vector<std::vector<mat4f>>& trajectories,
	std::vector<float>& baselines, unsigned int matchFileSkip = 1);

int main(int argc, char* argv[])
{
#if defined(DEBUG) | defined(_DEBUG)
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
	//_CrtSetBreakAlloc(1042);

	std::srand(0);

	try {
		//const std::string name = "matterport";	const std::string matchingPath = "//tirion/share/datasets/Matterport/Matching1/";	const std::string dataPath = "W:/data/matterport/v1_converted/";
		const std::string name = "scannet";	const std::string matchingPath = "//tirion/share/datasets/Matterport/MatchingScanNet/";	const std::string dataPath = "W:/data/scan-net/scans/checked/";
		const unsigned int matchFileSkip = 20;

		const std::string outFile = name + ".csv";
		const std::string trainTestPrefix = "scenes_";
		const std::vector<std::string> phases = { "test" }; //the current one is train/test for same scene

		std::ofstream ofs(outFile); const std::string splitter = ",";
		ofs << "scene" << splitter << "baseline values" << std::endl;
		for (const std::string& phase : phases) {
			const std::string fileList = matchingPath + trainTestPrefix + phase + ".txt";
			const std::vector<std::string> scenes = readScenes(fileList);
			//const std::vector<std::string> scenes = { "1LXtFkjw3qL" }; //debugging
			//const std::vector<std::string> scenes = { "2016-08-06_02-59-01__C7BA9586-8237-4204-9116-02AE5804338A" }; //debugging
			if (scenes.empty()) { std::cout << "warning: no scenes found from " << fileList << std::endl; continue; }
			std::unordered_map< std::string, std::vector<std::vector<mat4f>> > trajectories;
			std::vector<float> baselines;
			const std::string camAngleCacheFile = name + "-cameras.cache";
			if (util::fileExists(camAngleCacheFile)) {
				BinaryDataStreamFile s(camAngleCacheFile, false);
				s >> trajectories; s.close();
			}
			else {
				loadCameras(dataPath, scenes, trajectories);
				BinaryDataStreamFile s(camAngleCacheFile, true);
				s << trajectories; s.close();
			}
			float sumBaseline = 0.0f; unsigned int _idx = 0;
			std::cout << "computing.." << std::endl;
			for (const auto& a : trajectories) {
				std::vector<float> sceneBaselines;
				computeCameraBaselines(dataPath, matchingPath + a.first + "/matches1.txt", a.second, sceneBaselines, matchFileSkip);
				//write to file
				ofs << a.first << splitter;
				for (unsigned int i = 0; i < sceneBaselines.size(); i++) {
					ofs << sceneBaselines[i];
					if (i + 1 < sceneBaselines.size()) ofs << splitter;
				}
				ofs << std::endl;
				sumBaseline += std::accumulate(sceneBaselines.begin(), sceneBaselines.end(), 0.0f);
				baselines.insert(baselines.end(), sceneBaselines.begin(), sceneBaselines.end());
				std::cout << "\r[ " << ++_idx << " | " << trajectories.size() << " ]";
			}
			ofs.close();
			std::cout << std::endl;
			float average = sumBaseline / (float)baselines.size();
			std::cout << "average baseline: " << average << std::endl;
			std::vector<float> vdiff(baselines.size());
			std::transform(baselines.begin(), baselines.end(), vdiff.begin(), [average](float x) { return x - average; });
			float sqSum = std::inner_product(vdiff.begin(), vdiff.end(), vdiff.begin(), 0.0f);
			float stdev = std::sqrt(sqSum / (float)(baselines.size() - 1));
			std::cout << "stddev baseline: " << stdev << std::endl;
			std::cout << std::endl << "DONE DONE DONE" << std::endl;
		}
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

struct KeyMatch {
	vec2ui sensImIdx0;
	vec2ui sensImIdx1;
	vec2f pixPos0;
	vec2f pixPos1;
	vec3f worldPos0;
	vec3f worldPos1;
};
void readKeyMatches(const std::string& matchesFile, std::vector<KeyMatch>& res, unsigned int matchFileSkip)
{
	std::ifstream s(matchesFile);
	if (!s.is_open()) throw MLIB_EXCEPTION("failed to open file " + matchesFile + " for read");
	std::string line;
	size_t lineCount = 0, kpMatchCount = 0;
	while (std::getline(s, line)) {
		if (lineCount >= 3) { //skip header lines
			if ((lineCount - 3) % matchFileSkip == 0 || (lineCount - 3) % matchFileSkip == 1) {
				const auto parts = util::split(line, '\t');
				const unsigned int sensorIdx = util::convertTo<unsigned int>(parts[1]);
				const unsigned int imageIdx = util::convertTo<unsigned int>(parts[2]);
				vec2f pixelPos = util::convertTo<vec2f>(parts[3]);
				const vec3f worldPos = util::convertTo<vec3f>(parts[5]);

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

void computeCameraBaselines(const std::string& dataPath, const std::string& matchesFile, 
	const std::vector<std::vector<mat4f>>& trajectories,
	std::vector<float>& baselines, unsigned int matchFileSkip /*= 1*/)
{
	baselines.clear();
	std::vector<KeyMatch> keyMatches;
	readKeyMatches(matchesFile, keyMatches, matchFileSkip);
	std::cout << "#read poss = " << keyMatches.size() << std::endl;
	//compute baselines
	for (unsigned int i = 0; i < keyMatches.size(); i++) {
		const mat4f& cam0 = trajectories[keyMatches[i].sensImIdx0.x][keyMatches[i].sensImIdx0.y];
		const mat4f& cam1 = trajectories[keyMatches[i].sensImIdx1.x][keyMatches[i].sensImIdx1.y];
		float dist =  vec3f::dist(cam0.getTranslation(), cam1.getTranslation());
		if (dist > 0.1f)
			baselines.push_back(dist);
	}
}

void loadCameras(const std::string& dataPath, const std::vector<std::string>& scenes, std::unordered_map< std::string, std::vector<std::vector<mat4f>> >& trajectories)
{
	trajectories.clear();
	std::cout << "reading cameras from sens files..." << std::endl;

	unsigned int idx = 0;
	for (const std::string& scene : scenes) {
		const std::string scenePath = dataPath + scene + "/";
		const std::vector<std::string> sensFiles = Directory(scenePath).getFilesWithSuffix(".sens");
		if (sensFiles.size() > 3) throw MLIB_EXCEPTION("found " + std::to_string(sensFiles.size()) + " sens files for " + scene);
		std::vector<std::vector<mat4f>> trajectory(sensFiles.size());
		for (unsigned int i = 0; i < sensFiles.size(); i++) {
			SensorData sd(scenePath + sensFiles[i]);
			trajectory[i].resize(sd.m_frames.size());
			for (unsigned int f = 0; f < sd.m_frames.size(); f++)
				trajectory[i][f] = sd.m_frames[f].getCameraToWorld();
		}
		trajectories[scene] = trajectory;
		std::cout << "\r[ " << ++idx << " | " << scenes.size() << " ]";
	}
	std::cout << "done!" << std::endl;
}


