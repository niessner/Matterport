
#pragma once

#include "mLibInclude.h"

class NormalExtractor {
public:

	NormalExtractor(unsigned int width, unsigned int height) {
		m_width = width;
		m_height = height;
		m_graphics = new D3D11GraphicsDevice();
		m_graphics->initWithoutWindow();

		m_cbCamera.init(*m_graphics);
		m_shaders.init(*m_graphics);
		m_shaders.registerShader("shaders/normals.hlsl", "normals");

		std::vector<DXGI_FORMAT> formats = {
			DXGI_FORMAT::DXGI_FORMAT_R32G32B32A32_FLOAT
		};
		m_renderTarget.init(*m_graphics, m_width, m_height, formats);
	}
	~NormalExtractor() {
		SAFE_DELETE(m_graphics);
	}

	//normals from bilateral filt depth
	void computeDepthNormals(const SensorData* sd, float depthSigmaD, float depthSigmaR, std::vector<PointImage>& normals);

	void computeMeshNormals(const SensorData* sd, const TriMeshf& mesh, float zNear, float zFar, std::vector<PointImage>& normals);


	static void saveNormalImage(const std::string& filename, const PointImage& normal)
	{
		//convert to 3x 16bit
		BaseImage<vec3us> normalImage(normal.getDimensions());
		for (auto& p : normal) {
			vec3f n = (p.value.x == -std::numeric_limits<float>::infinity()) ? vec3f(0.0f, 0.0f, 0.0f) : p.value;
			n = 0.5f * (n + 1) * 65535.0f; //[-1,1] -> [0,65535]
			normalImage(p.x, p.y) = vec3us(n);
		}

		FreeImageWrapper::saveImage(filename, normalImage);
	}
	static void loadNormalImage(const std::string& filename, PointImage& normal)
	{
		//convert to 3x 16bit
		BaseImage<vec3us> normalImage;
		FreeImageWrapper::loadImage(filename, normalImage);
		normal.allocate(normalImage.getWidth(), normalImage.getHeight());

		for (auto& p : normal) {
			vec3f n = vec3f(normalImage(p.x, p.y)) * 2.0f / 65535.0f - 1.0f;
			if (n.length() < 0.001f)
				n = vec3f(-std::numeric_limits<float>::infinity());
			p.value = n;
		}
	}


private:
	static mat4f projFromVision(const mat4f& intrinsic, const mat4f& cameraToWorld, unsigned int width, unsigned int height, float zNear, float zFar);
	void render(const mat4f& perspective, const mat4f& worldToCam, const vec3f& eye, ColorImageR32G32B32A32& image, const SensorData* sd)
	{
		m_renderTarget.bind();
		m_renderTarget.clear(ml::vec4f(0.0f, 0.0f, 0.0f, 0.0f));

		ConstantBufferCamera cbCamera;
		cbCamera.worldViewProj = perspective * worldToCam;
		cbCamera.world = worldToCam;
		m_cbCamera.updateAndBind(cbCamera, 0);
		m_shaders.bindShaders("normals");
		m_mesh.render();

		m_renderTarget.unbind();
		m_renderTarget.captureColorBuffer(image);

		//mat4f intrinsic = sd->m_calibrationDepth.m_intrinsic;
		//if (m_width != sd->m_depthWidth || m_height != sd->m_depthHeight) {
		//	// adapt depth intrinsics
		//	intrinsic._m00 *= (float)m_width / (float)sd->m_depthWidth;				//focal length
		//	intrinsic._m11 *= (float)m_height / (float)sd->m_depthHeight;			//focal length
		//	intrinsic._m02 *= (float)(m_width - 1) / (float)(sd->m_depthWidth - 1);		//principal point
		//	intrinsic._m12 *= (float)(m_height - 1) / (float)(sd->m_depthHeight - 1);	//principal point
		//}
		//DepthImage32 depth;
		//m_renderTarget.captureDepthBuffer(depth);
		//mat4f projToCamera = perspective.getInverse();
		//for (auto &p : depth) {
		//	vec3f posWorld = vec3f(-std::numeric_limits<float>::infinity());
		//	if (p.value != 0.0f && p.value != 1.0f) {
		//		vec3f posProj = vec3f(m_graphics->pixelToNDC(vec2i((int)p.x, (int)p.y), depth.getWidth(), depth.getHeight()), p.value);
		//		vec3f posCamera = projToCamera * posProj;

		//		if (posCamera.z >= 0.4f && posCamera.z <= 7.0f) {
		//			p.value = posCamera.z;
		//			posWorld = worldToCam.getInverse() * posCamera;
		//		}
		//		else {
		//			p.value = -std::numeric_limits<float>::infinity();
		//		}
		//	}
		//	else {
		//		p.value = -std::numeric_limits<float>::infinity();
		//	}
		//} //depth pixels
		//PointImage normal = sd->computeNormals(intrinsic.getInverse(), depth);
		//PointImage campos = sd->computeCameraSpacePositions(0);
		//campos.resize(m_width, m_height);
		//{
		//	PointCloudf pc;
		//	for (const auto& p : campos) {
		//		vec3f n = normal(p.x, p.y);
		//		if (p.value.x != -std::numeric_limits<float>::infinity() && (!(n.x == 0 && n.y == 0 && n.z == 0))) {
		//			pc.m_points.push_back(p.value);
		//			pc.m_normals.push_back(n);
		//		}
		//	}
		//	PointCloudIOf::saveToFile("test-d.ply", pc);
		//}
		//saveNormalImage("test-d.png", normal);
		//int a = 5;
	}

	struct ConstantBufferCamera {
		mat4f worldViewProj;
		mat4f world;
		//vec4f eye;
	};

	unsigned int m_width;
	unsigned int m_height;

	D3D11GraphicsDevice* m_graphics;
	D3D11ShaderManager m_shaders;
	D3D11ConstantBuffer<ConstantBufferCamera> m_cbCamera;
	D3D11TriMesh m_mesh;
	D3D11RenderTarget m_renderTarget;
};
