// sensorFile.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iomanip>


template <typename T> BaseImage<T> undistort(const BaseImage<T>& src, const mat4f& intrinsic, const float coeff[5])
{
	BaseImage<T> res(src.getWidth(), src.getHeight());
	res.setInvalidValue(src.getInvalidValue());
	res.setPixels(res.getInvalidValue());

	for (unsigned int y = 0; y < src.getHeight(); y++)	{
		for (unsigned int x = 0; x < src.getWidth(); x++)	{
			vec2f nic_loc;
			vec2f sample_loc;

			//Normalized image coords
			nic_loc.x = (x - intrinsic(0, 2)) / intrinsic(0, 0);
			nic_loc.y = (y - intrinsic(1, 2)) / intrinsic(1, 1);

			float r2 = nic_loc.x * nic_loc.x + nic_loc.y * nic_loc.y;

			// Radial distortion
			sample_loc.x = nic_loc.x * (1.0f + r2 * coeff[0] + r2*r2 * coeff[1] + r2*r2*r2 * coeff[4]);
			sample_loc.y = nic_loc.y * (1.0f + r2 * coeff[0] + r2*r2 * coeff[1] + r2*r2*r2 * coeff[4]);

			// Tangential distortion
			sample_loc.x += 2.0f * coeff[2] * nic_loc.x * nic_loc.y + coeff[3] * (r2 + 2.0f * nic_loc.x * nic_loc.x);
			sample_loc.y += coeff[2] * (r2 + 2.0f * nic_loc.y * nic_loc.y) + 2.0f * coeff[3] * nic_loc.x * nic_loc.y;

			// Move back to the image space
			sample_loc.x = sample_loc.x * intrinsic(0, 0) + intrinsic(0, 2);
			sample_loc.y = sample_loc.y * intrinsic(1, 1) + intrinsic(1, 2);

			vec2i sample_loc_i = math::round(sample_loc);
			if (src.isValidCoordinate(sample_loc_i)) {
				res(x, y) = src(sample_loc_i);
			}
		}
	}
	return res;
}

void convertToSens(const std::string& srcPath, const std::string& outFile) 
{
	Directory dir(srcPath);
	std::vector<std::string> baseFiles = dir.getFilesWithSuffix("_d0_0.png");
	for (auto& f : baseFiles) f = util::replace(f, "_d0_0.png", "");

	const unsigned int numCams = 3;
	SensorData sd[numCams];
	for (size_t fidx = 0; fidx < baseFiles.size(); fidx++) {
		std::cout << "\rProcessed [ " << fidx << " | " << baseFiles.size() << " ] ";
		const std::string& f = baseFiles[fidx];		
		for (unsigned int camIdx = 0; camIdx < numCams; camIdx++) {
			for (unsigned int poseIdx = 0; poseIdx < 6; poseIdx++) {
				const std::string depthFile = srcPath + "/" + f + "_d" + std::to_string(camIdx) + "_" + std::to_string(poseIdx) + ".png";
				const std::string colorFile = srcPath + "/" + f + "_i" + std::to_string(camIdx) + "_" + std::to_string(poseIdx) + ".jpg";
				const std::string poseFile = srcPath + "/" + f + "_pose_" + std::to_string(camIdx) + "_" + std::to_string(poseIdx) + ".txt";
				const std::string intrFile = srcPath + "/" + f + "_intrinsics_" + std::to_string(camIdx) + ".txt";

				DepthImage16 depthImage16;
				FreeImageWrapper::loadImage(depthFile, depthImage16);
				depthImage16.setInvalidValue(0);
				//for (auto& p : depthImage16) p.value /= 5;

				//DepthImage32 depthImage32(depthImage16.getWidth(), depthImage16.getHeight());
				//depthImage32.setInvalidValue(0.0f);
				//for (auto& p : depthImage32) p.value = 0.001f * depthImage16(p.x, p.y);
				//FreeImageWrapper::saveImage("test.png", ColorImageR32G32B32A32(depthImage32, true));
				//std::cout << "DONE" << std::endl;
				//getchar();
				//exit(1);

				ColorImageR8G8B8 colorImage;
				FreeImageWrapper::loadImage(colorFile, colorImage);

				mat4f pose;
				{
					std::ifstream inFile(poseFile);
					for (unsigned int i = 0; i < 16; i++) inFile >> pose[i];
				}
				unsigned int width, height;
				float fx, fy, mx, my;
				float k[5];
				{
					std::ifstream inFile(intrFile);
					inFile >> width >> height;
					inFile >> fx >> fy >> mx >> my;
					for (unsigned i = 0; i < 5; i++) inFile >> k[i];
				}

				if (fidx == 0 && poseIdx == 0) {
					SensorData::CalibrationData calibrationDepth(
						mat4f(
						fx, 0.0f, mx, 0.0f,
						0.0f, fy, my, 0.0f,
						0.0f, 0.0f, 1.0f, 0.0f,
						0.0f, 0.0f, 0.0f, 1.0f
						));
					SensorData::CalibrationData calibrationColor = calibrationDepth;

					float depthShift = 1.0f / 0.00025f;	//depth values are 0.25mm
					sd[camIdx].initDefault(width, height, width, height, calibrationColor, calibrationDepth,
						SensorData::COMPRESSION_TYPE_COLOR::TYPE_JPEG, SensorData::COMPRESSION_TYPE_DEPTH::TYPE_ZLIB_USHORT, depthShift, "Matterport");
				}

				colorImage = undistort(colorImage, sd[camIdx].m_calibrationColor.m_intrinsic, k);
				depthImage16 = DepthImage16(undistort(depthImage16, sd[camIdx].m_calibrationDepth.m_intrinsic, k));
				sd[camIdx].addFrame(colorImage.getData(), depthImage16.getData(), pose);
			}
		}
	}
	std::cout << std::endl;
	
	for (unsigned int i = 0; i < numCams; i++) {
		//std::cout << sd[i] << std::endl;

		//sd[i].saveToFile("test_" + std::to_string(i) + ".sens");
		auto s = util::splitOnLast(outFile, ".");
		std::string _outFile = s.first + "_" + std::to_string(i) + "." + s.second;
		sd[i].saveToFile(_outFile);
	}
}

int main(int argc, char* argv[])
{
#if defined(DEBUG) | defined(_DEBUG)
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
	//_CrtSetBreakAlloc(7545);
	try {
		const std::string srcFolder = "W:/data/matterport/v1";
		const std::string dstFolder = "W:/data/matterport/v1_converted";

		util::makeDirectory(dstFolder);

		Directory srcDir(srcFolder);

		for (const std::string& s : srcDir.getDirectories()) {
			std::cout << s << std::endl;			
			
			const std::string outPath = dstFolder + "/" + s;
			if (true) {
				util::deleteDirectory(outPath);
				util::makeDirectory(outPath);
				if (true) {
					//extract the mesh
					const std::string srcFileMesh = srcFolder + "/" + s + "/" + s + "_mesh.zip";
					const std::string cmd = "unzip -qq " + srcFileMesh + " -d " + outPath;
					system(cmd.c_str());
					const std::vector<std::string> dirs = Directory::enumerateDirectories(outPath);
					if (dirs.size() != 1) throw MLIB_EXCEPTION("excepting exactly one output directory");
					const std::string subDirSrc = outPath + "/" + dirs.front();
					const std::string newSubDirSrc = outPath + "/mesh/";
					util::renameFile(subDirSrc, newSubDirSrc);
				}
				if (true) {
					//extract the raw data
					const std::string srcFileData = srcFolder + "/" + s + "/" + s + ".zip";
					const std::string cmd = "unzip -qq " + srcFileData + " -d " + outPath;
					system(cmd.c_str());
					const std::string subDirSrc = outPath + "/" + s;
					const std::string newSubDirSrc = outPath + "/data/";
					util::renameFile(subDirSrc, newSubDirSrc);
				}
			}

			convertToSens(outPath + "/data/", outPath + "/" + s + ".sens");
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

