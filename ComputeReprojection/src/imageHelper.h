#pragma once

class ImageHelper
{
public:
	ImageHelper() {}
	~ImageHelper() {}


	static inline float gaussR(float sigma, float dist)
	{
		return exp(-(dist*dist) / (2.0f*sigma*sigma));
	}
	static inline float gaussR(float sigma, const vec3f& d)
	{
		float dist = d.length();
		return exp(-(dist*dist) / (2.0f*sigma*sigma));
	}
	static inline float gaussR(float sigma, const vec3uc& d)
	{
		vec3f _d(d);	//_d /= 255.0f;
		float dist = _d.length();
		return exp(-(dist*dist) / (2.0f*sigma*sigma));
	}


	static inline float linearR(float sigma, float dist)
	{
		return std::max(1.0f, std::min(0.0f, 1.0f - (dist*dist) / (2.0f*sigma*sigma)));
	}

	static inline float gaussD(float sigma, int x, int y)
	{
		return exp(-((x*x + y*y) / (2.0f*sigma*sigma)));
	}

	static inline float gaussD(float sigma, int x)
	{
		return exp(-((x*x) / (2.0f*sigma*sigma)));
	}

	static void gaussFilter(BaseImage<float>& img, float sigmaD) {

		BaseImage<float> res(img.getDimensions());
		res.setInvalidValue(img.getInvalidValue());

		const int kernelRadius = (int)ceil(2.0*sigmaD);
#pragma omp parallel for
		for (int _y = 0; _y < (int)img.getHeight(); _y++) {
			unsigned int y = (unsigned int)_y;
			for (unsigned int x = 0; x < img.getWidth(); x++) {

				res.setInvalid(x, y);

				float sum = 0.0f;
				float sumWeight = 0.0f;

				if (img.isValid(x, y)) {
					const float center = img(x, y);

					for (int m = x - kernelRadius; m <= (int)x + kernelRadius; m++) {
						for (int n = y - kernelRadius; n <= (int)y + kernelRadius; n++) {
							if (m >= 0 && n >= 0 && m < (int)img.getWidth() && n < (int)img.getHeight()) {
								if (img.isValid(m, n)) {
									const float current = img(m, n);
									const float weight = gaussD(sigmaD, m - x, n - y);
									sumWeight += weight;
									sum += weight*current;
								}
							}
						}
					}
					if (sumWeight > 0.0f) res(x, y) = sum / sumWeight;
				}
			}
		}
		img = res;
	}
	static void bilateralFilter(DepthImage32& d, float sigmaD, float sigmaR) {

		DepthImage32 res(d.getWidth(), d.getHeight());
		res.setInvalidValue(d.getInvalidValue());

		const int kernelRadius = (int)ceil(2.0*sigmaD);
#pragma omp parallel for
		for (int _y = 0; _y < (int)d.getHeight(); _y++) {
			unsigned int y = (unsigned int)_y;
			for (unsigned int x = 0; x < d.getWidth(); x++) {
				res.setInvalid(x, y);

				float sum = 0.0f;
				float sumWeight = 0.0f;

				if (d.isValid(x, y)) {
					const float center = d(x, y);

					for (int m = x - kernelRadius; m <= (int)x + kernelRadius; m++) {
						for (int n = y - kernelRadius; n <= (int)y + kernelRadius; n++) {
							if (m >= 0 && n >= 0 && m < (int)d.getWidth() && n < (int)d.getHeight()) {
								if (d.isValid(m, n)) {
									const float current = d(m, n);
									const float weight = gaussD(sigmaD, m - x, n - y)*gaussR(sigmaR, current - center);
									sumWeight += weight;
									sum += weight*current;
								}
							}
						}
					}
					if (sumWeight > 0.0f) res(x, y) = sum / sumWeight;
				}
			}
		}
		d = res;
	}

	static void erode(DepthImage32& depth, unsigned int numIter = 2) {
		numIter = 2 * ((numIter + 1) / 2);
		DepthImage32 tmp; tmp.setInvalidValue(depth.getInvalidValue());
		for (unsigned int i = 0; i < numIter; i++) {
			if (i % 2 == 0) {
				erode(tmp, depth, 3, 0.05f, 0.3f);
			}
			else {
				erode(depth, tmp, 3, 0.05f, 0.3f);
			}
		}
	}

	static BaseImage<float> convertToGrayscale(const ColorImageR8G8B8& image) {
		BaseImage<float> res(image.getWidth(), image.getHeight());
		res.setInvalidValue(0.0f);

		for (const auto& p : image) {
			float v = (0.299f*p.value.x + 0.587f*p.value.y + 0.114f*p.value.z) / 255.0f;
			res(p.x, p.y) = v;
		}

		return res;
	}


	static BaseImage<float> computeGradientMagnitude(const ColorImageR32& image)
	{
		BaseImage<float> res(image.getWidth(), image.getHeight());
		res.setInvalidValue(-std::numeric_limits<float>::infinity());
		res.setPixels(res.getInvalidValue());
		const auto invalid = image.getInvalidValue();

#pragma omp parallel for
		for (int _y = 0; _y < (int)image.getHeight(); _y++) {
			unsigned int y = (unsigned int)_y;
			for (unsigned int x = 0; x < image.getWidth(); x++) {
				if (x > 0 && x < image.getWidth() - 1 && y > 0 && y < image.getHeight() - 1) {
					float pos00 = image(x - 1, y - 1);	if (pos00 == invalid) continue;
					float pos01 = image(x - 1, y - 0);	if (pos01 == invalid) continue;
					float pos02 = image(x - 1, y + 1);	if (pos02 == invalid) continue;

					float pos10 = image(x - 0, y - 1); if (pos10 == invalid) continue;
					float pos12 = image(x - 0, y + 1); if (pos12 == invalid) continue;

					float pos20 = image(x + 1, y - 1); if (pos20 == invalid) continue;
					float pos21 = image(x + 1, y - 0); if (pos21 == invalid) continue;
					float pos22 = image(x + 1, y + 1); if (pos22 == invalid) continue;

					float resU = (-1.0f)*pos00 + (1.0f)*pos20 +
						(-2.0f)*pos01 + (2.0f)*pos21 +
						(-1.0f)*pos02 + (1.0f)*pos22;
					resU /= 8.0f;

					float resV = (-1.0f)*pos00 + (-2.0f)*pos10 + (-1.0f)*pos20 +
						(1.0f)*pos02 + (2.0f)*pos12 + (1.0f)*pos22;
					resV /= 8.0f;

					res(x, y) = vec2f(resU, resV).length();
				}
			}
		}
		return res;
	}

	static void computeImageStatistics(const BaseImage<float>& image)
	{
		float min = std::numeric_limits<float>::infinity(), max = -std::numeric_limits<float>::infinity();
		float mean = 0.0f; unsigned int count = 0;
		for (const auto& p : image) {
			if (p.value != image.getInvalidValue()) {
				if (p.value < min) min = p.value;
				if (p.value > max) max = p.value;
				mean += p.value;
				count++;
			}
		}
		mean /= count;
		std::cout << "image range [" << min << ", " << max << "]" << std::endl;
		std::cout << "mean value = " << mean << " (" << count << "/" << image.getNumPixels() << " valid pixels)" << std::endl;
	}

private:

	static void erode(DepthImage32& output, const DepthImage32& input, int structureSize, float dThresh, float fracReq)
	{
		output.allocate(input.getWidth(), input.getHeight());

		for (unsigned int y = 0; y < input.getHeight(); y++) {
			for (unsigned int x = 0; x < input.getWidth(); x++) {
				unsigned int count = 0;
				float oldDepth = input(x, y);
				for (int i = -structureSize; i <= structureSize; i++) {
					for (int j = -structureSize; j <= structureSize; j++) {
						if (x + j >= 0 && x + j < input.getWidth() && y + i >= 0 && y + i < input.getHeight()) {
							float depth = input(x + j, y + i);
							if (depth == input.getInvalidValue() || depth == 0.0f || fabs(depth - oldDepth) > dThresh) {
								count++;
							}
						}
					}
				}

				unsigned int sum = (2 * structureSize + 1)*(2 * structureSize + 1);
				if ((float)count / (float)sum >= fracReq) {
					output(x, y) = input.getInvalidValue();
				}
				else {
					output(x, y) = input(x, y);
				}
			}
		}
	}
};

