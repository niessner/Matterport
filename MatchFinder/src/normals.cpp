
#include "stdafx.h"

#include "normals.h"
#include "imageHelper.h"


void NormalExtractor::computeDepthNormals(const SensorData* sd, const float depthSigmaD, const float depthSigmaR,
	unsigned int outWidth, unsigned int outHeight, std::vector<PointImage>& normals)
{
	const mat4f depthIntrinsicInv = sd->m_calibrationDepth.m_intrinsic.getInverse();
	normals.resize(sd->m_frames.size());

	for (size_t imageIdx = 0; imageIdx < sd->m_frames.size(); imageIdx++) {
		DepthImage32 d = sd->computeDepthImage(imageIdx);
		d.resize(outWidth, outHeight);
		DepthImage32 dfilt = d; 
		ImageHelper::bilateralFilter(dfilt, depthSigmaD, depthSigmaR);
		normals[imageIdx] = SensorData::computeNormals(depthIntrinsicInv, dfilt);
	}
}
