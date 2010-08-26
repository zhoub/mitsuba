#include <mitsuba/core/plugin.h>
#include <mitsuba/hw/vpl.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/hw/gputexture.h>

MTS_NAMESPACE_BEGIN

VPLShaderManager::VPLShaderManager(const Scene *scene, Renderer *renderer)
	 : m_scene(scene), m_renderer(renderer), m_clamping(0.1f),
 	   m_maxClipDist(std::numeric_limits<Float>::infinity()), m_initialized(false), 
	   m_shadowMapResolution(512), m_singlePass(false), m_allowNonDiffuseVPLs(false) {
}

VPLShaderManager::~VPLShaderManager() {
	if (m_initialized)
		cleanup();
}

void VPLShaderManager::init() {
    if (m_renderer->getCapabilities()->isSupported(RendererCapabilities::EGeometryShaders)) {
		m_shadowProgram = m_renderer->createGPUProgram("Shadow Program");
	
		m_shadowProgram->setSource(GPUProgram::EVertexProgram,
			"void main() {\n"
			"	gl_Position = gl_Vertex;\n"
			"}"
		);

		m_shadowProgram->setSource(GPUProgram::EGeometryProgram,
			"#version 120\n"
			"#extension GL_EXT_geometry_shader4 : enable\n"
			"\n"
			"uniform mat4 cubeMapTransform[6];\n"
			"uniform vec4 depthVec[6];\n"
			"varying float depth;\n"
			"\n"
			"void main() {\n"
			"	depth = 0;\n" // avoid a (incorrect?) warning
			"   for (int side = 0; side < 6; side++) {\n"
			"	    gl_Layer = side;\n"
			"       for (int i = 0; i < gl_VerticesIn; i++) {\n"
			"           gl_Position = cubeMapTransform[side] * gl_PositionIn[i];\n"
			"           depth = dot(depthVec[side], gl_PositionIn[i]);\n"
			"           EmitVertex();\n"
			"       }\n"
			"       EndPrimitive();\n"
			"   }\n"
			"}\n"
		);

		m_shadowProgram->setSource(GPUProgram::EFragmentProgram,
			"#version 120\n"
			"varying float depth;\n"
			"void main() {\n"
			"	gl_FragDepth = depth;\n"
			"	float dx = dFdx(depth), dy = dFdy(depth);"
			"   gl_FragDepth = depth + sqrt(dx*dx + dy*dy);"
			"}\n"
		);
	
		/* Six output triangles per input triangle */
		m_shadowProgram->setMaxVertices(18); 
		m_shadowProgram->init();

		for (int i=0; i<6; ++i)	{
			m_shadowProgramParam_cubeMapTransform[i] = 
				m_shadowProgram->getParameterID(formatString("cubeMapTransform[%i]", i));
			m_shadowProgramParam_depthVec[i] = 
				m_shadowProgram->getParameterID(formatString("depthVec[%i]", i));
		}
	}

	m_altShadowProgram = m_renderer->createGPUProgram("Alternative Shadow Program");
	m_altShadowProgram->setSource(GPUProgram::EVertexProgram,
		"uniform mat4 cubeMapTransform;\n"
		"uniform vec4 depthVec;\n"
		"varying float depth;\n"
		"void main() {\n"
		"	gl_Position = cubeMapTransform * gl_Vertex;\n"
		"	depth = dot(depthVec, gl_Vertex);\n"
		"}\n"
	);

	m_altShadowProgram->setSource(GPUProgram::EFragmentProgram,
		"#version 120\n"
		"varying float depth;\n"
		"void main() {\n"
		"	float dx = dFdx(depth), dy = dFdy(depth);"
		"   gl_FragDepth = depth + sqrt(dx*dx + dy*dy);"
		"}\n"
	);

	m_altShadowProgram->init();

	m_altShadowProgramParam_cubeMapTransform = 
		m_altShadowProgram->getParameterID("cubeMapTransform");
	m_altShadowProgramParam_depthVec = 
		m_altShadowProgram->getParameterID("depthVec");

	const std::vector<TriMesh *> meshes = m_scene->getMeshes();
	const std::vector<Luminaire *> luminaires = m_scene->getLuminaires();

	for (size_t i=0; i<meshes.size(); ++i) {
		m_renderer->registerGeometry(meshes[i]);
		Shader *shader = m_renderer->registerShaderForResource(meshes[i]->getBSDF());
		if (shader != NULL && !shader->isComplete())
			m_renderer->unregisterShaderForResource(meshes[i]->getBSDF());
	}

	for (size_t i=0; i<luminaires.size(); ++i)
		m_renderer->registerShaderForResource(luminaires[i]);

	if (m_scene->hasBackgroundLuminaire() && 
		m_renderer->getShaderForResource(m_scene->getBackgroundLuminaire()) != NULL) {
		Shader *shader = m_renderer->getShaderForResource(m_scene->getBackgroundLuminaire());
		m_backgroundDependencies = VPLDependencyNode(shader);
		int id = 0;
		std::ostringstream oss;
		std::string evalName = m_backgroundDependencies.recursiveGenerateCode(oss, id);

		m_backgroundProgram = m_renderer->createGPUProgram("Background program");
		m_backgroundProgram->setSource(GPUProgram::EVertexProgram,
			"uniform mat4 clipToWorld;\n"
			"varying vec3 d;\n"
			"void main() {\n"
			"	gl_Position = ftransform();\n"
			"	vec4 tmp = clipToWorld * (gl_ModelViewProjectionMatrix * gl_Vertex);\n"
			"   d = tmp.xyz/tmp.w;"
			"}\n"
		);

		oss << "varying vec3 d;" << endl
			<< "uniform vec3 camPos;" << endl
			<< "void main() {" << endl
			<< "  gl_FragColor.rgb = " << evalName << "_background(normalize(d - camPos));" << endl
			<< "  gl_FragColor.a = 1.0;" << endl
			<< "}" << endl;

		m_backgroundProgram->setSource(GPUProgram::EFragmentProgram, oss.str());
		m_backgroundProgram->init();

		id = 0;
		m_backgroundDependencies.recursiveResolve(m_backgroundProgram, id);
	}

	m_initialized = true;
}

void VPLShaderManager::setVPL(const VPL &vpl) {
	Point p = vpl.its.p + vpl.its.shFrame.n * 0.01;
	Intersection its;

	/* Estimate good near and far plane locations by tracing some rays */
	Float nearClip =  std::numeric_limits<Float>::infinity(),
		  farClip  = -std::numeric_limits<Float>::infinity();
	Ray ray;
	ray.o = p;

	if (m_shadowMap == NULL || m_shadowMapResolution != m_shadowMap->getSize().x) {
		m_shadowMap = m_renderer->createGPUTexture("Shadow cube map", NULL);
		m_shadowMap->setSize(Point3i(m_shadowMapResolution, m_shadowMapResolution, 1));
		m_shadowMap->setFrameBufferType(GPUTexture::EDepthBuffer);
		m_shadowMap->setType(GPUTexture::ETextureCubeMap);
		m_shadowMap->setWrapType(GPUTexture::EClampToEdge);
		m_shadowMap->setFilterType(GPUTexture::ENearest);
		m_shadowMap->setDepthMode(GPUTexture::ENormal);
		m_shadowMap->init();
	}

	const int sampleCount = 200;
	const Float invSampleCount = 1.0f/sampleCount;

	for (int i=1; i<=sampleCount; ++i) {
		Vector dir;
		Point2 seed(i*invSampleCount, radicalInverse(2, i)); // Hammersley seq.
		if (vpl.type == ELuminaireVPL && vpl.luminaire->getType() & Luminaire::EOnSurface)
			dir = squareToSphere(seed);
		else
			dir = vpl.its.shFrame.toWorld(squareToHemispherePSA(seed));
		ray.setDirection(dir);
		if (m_scene->rayIntersect(ray, its)) {
			nearClip = std::min(nearClip, its.t);
			farClip = std::max(farClip, its.t);
		}
	}

	m_minDist = nearClip + (farClip - nearClip) * m_clamping;

	nearClip = std::min(nearClip, (Float) 0.001f);
	farClip = std::min(farClip * 1.5f, m_maxClipDist);

	if (farClip < 0 || nearClip >= farClip) {
		/* Unable to find any surface - just default values based on the scene size */
		nearClip = 1e-3 * m_scene->getBSphere().radius;
		farClip = 1e3 * m_scene->getBSphere().radius;
		m_minDist = 0;
	}

	m_nearClip = nearClip;
	m_invClipRange = 1/(farClip-nearClip);
	Transform lightViewTrafo, lightProjTrafo = Transform::glPerspective(90.0f, nearClip, farClip);

	m_shadowMap->activateTarget();
	if (m_singlePass && m_shadowProgram != NULL) {
		/* "Fancy": render the whole cube map in a single pass using 
		   a geometry program. On anything but brand-new hardware, this 
		   is actually slower. */

		m_shadowMap->activateSide(-1);
		m_shadowMap->clear();
		m_shadowProgram->bind();
		for (int i=0; i<6; ++i) {
			switch (i) {
				case 0: lightViewTrafo = Transform::lookAt(p, p + Vector(1, 0, 0), Vector(0, 1, 0)).inverse(); break;
				case 1: lightViewTrafo = Transform::lookAt(p, p + Vector(-1, 0, 0), Vector(0, 1, 0)).inverse(); break;
				case 2: lightViewTrafo = Transform::lookAt(p, p + Vector(0, 1, 0), Vector(0, 0, -1)).inverse(); break;
				case 3: lightViewTrafo = Transform::lookAt(p, p + Vector(0, -1, 0), Vector(0, 0, 1)).inverse(); break;
				case 4: lightViewTrafo = Transform::lookAt(p, p + Vector(0, 0, 1), Vector(0, 1, 0)).inverse(); break;
				case 5: lightViewTrafo = Transform::lookAt(p, p + Vector(0, 0, -1), Vector(0, 1, 0)).inverse(); break;
			}
			const Matrix4x4 *viewMatrix = lightViewTrafo.getMatrix();
			m_shadowProgram->setParameter(m_shadowProgramParam_cubeMapTransform[i], lightProjTrafo * lightViewTrafo);
			m_shadowProgram->setParameter(m_shadowProgramParam_depthVec[i], Vector4(
				-viewMatrix->m[2][0] * m_invClipRange,
				-viewMatrix->m[2][1] * m_invClipRange,
				-viewMatrix->m[2][2] * m_invClipRange,
				(-viewMatrix->m[2][3] - m_nearClip) * m_invClipRange
			));
		}
		m_renderer->drawAll();
		m_shadowProgram->unbind();
	} else {
		/* Old-fashioned: render 6 times, once for each cube map face */

		m_altShadowProgram->bind();
		for (int i=0; i<6; ++i) {
			switch (i) {
				case 0: lightViewTrafo = Transform::lookAt(p, p + Vector(1, 0, 0), Vector(0, 1, 0)).inverse(); break;
				case 1: lightViewTrafo = Transform::lookAt(p, p + Vector(-1, 0, 0), Vector(0, 1, 0)).inverse(); break;
				case 2: lightViewTrafo = Transform::lookAt(p, p + Vector(0, 1, 0), Vector(0, 0, -1)).inverse(); break;
				case 3: lightViewTrafo = Transform::lookAt(p, p + Vector(0, -1, 0), Vector(0, 0, 1)).inverse(); break;
				case 4: lightViewTrafo = Transform::lookAt(p, p + Vector(0, 0, 1), Vector(0, 1, 0)).inverse(); break;
				case 5: lightViewTrafo = Transform::lookAt(p, p + Vector(0, 0, -1), Vector(0, 1, 0)).inverse(); break;
			}
			const Matrix4x4 *viewMatrix = lightViewTrafo.getMatrix();

			m_altShadowProgram->setParameter(m_altShadowProgramParam_cubeMapTransform, lightProjTrafo * lightViewTrafo);
			m_altShadowProgram->setParameter(m_altShadowProgramParam_depthVec, Vector4(
				-viewMatrix->m[2][0] * m_invClipRange,
				-viewMatrix->m[2][1] * m_invClipRange,
				-viewMatrix->m[2][2] * m_invClipRange,
				(-viewMatrix->m[2][3] - m_nearClip) * m_invClipRange
			));

			m_shadowMap->activateSide(i);
			m_shadowMap->clear();
			m_renderer->drawAll();
		}

		m_altShadowProgram->unbind();
	}
	m_shadowMap->releaseTarget();
}

void VPLShaderManager::configure(const VPL &vpl, const BSDF *bsdf, const Luminaire *luminaire, const Point &camPos) {
	Shader *bsdfShader = m_renderer->getShaderForResource(bsdf);
	Shader *vplShader = (vpl.type == ELuminaireVPL)
		? m_renderer->getShaderForResource(vpl.luminaire)
		: m_renderer->getShaderForResource(vpl.its.shape->getBSDF());
	Shader *lumShader = (luminaire == NULL) ? NULL 
		: m_renderer->getShaderForResource(luminaire);
	std::ostringstream oss;

	
	if (bsdfShader == NULL || vplShader == NULL ||
		(luminaire != NULL && lumShader == NULL)) {
		/* Unsupported! */
		m_renderer->setColor(Spectrum(0.0f));
		return;
	}

	m_targetConfig = VPLProgramConfiguration(vplShader, bsdfShader, lumShader);
	m_targetConfig.toString(oss);
	std::string configName = oss.str();
	std::map<std::string, ProgramAndConfiguration>::iterator it =
		m_programs.find(configName);
	GPUProgram *program = NULL;

	if (it != m_programs.end()) {
		/* A program for this configuration has been created previously */
		m_current = (*it).second;
		program = m_current.program;
	} else {
		/* No program for this particular combination exists -- create one */
		program = m_renderer->createGPUProgram(configName);

		/* Vertex program */
		oss.str("");
        oss << "#version 120" << endl
			<< "uniform vec3 vplPos, camPos;" << endl
			<< "varying vec3 normal, tangent, lightVec, camVec;" << endl
			<< "varying vec2 uv;" << endl
			<< endl
			<< "void main() {" << endl
			<< "   normal = gl_Normal;" << endl
			<< "   tangent = gl_MultiTexCoord1.xyz;" << endl
			<< "   uv = gl_MultiTexCoord0.xy;" << endl
			<< "   camVec = camPos - gl_Vertex.xyz;" << endl
			<< "   lightVec = vplPos - gl_Vertex.xyz;" << endl
            << "   gl_Position = ftransform();" << endl
			<< "}" << endl;

		program->setSource(GPUProgram::EVertexProgram, oss.str());
		oss.str("");

		oss << "#version 120" << endl
			<< endl
			<< "/* Uniform inputs */" << endl
			<< "uniform samplerCube shadowMap;" << endl
			<< "uniform vec3 vplPower, vplS, vplT, vplN, vplWi;" << endl
			<< "uniform float nearClip, invClipRange, minDist;" << endl
			<< "uniform vec2 vplUV;" << endl
			<< "uniform bool allowNonDiffuseVPLs;" << endl
			<< endl
			<< "/* Inputs <- Vertex program */" << endl
			<< "varying vec3 normal, tangent, lightVec, camVec;" << endl
			<< "varying vec2 uv;" << endl
			<< endl
			<< "/* Some helper functions for BSDF implementations */" << endl
			<< "float cosTheta(vec3 v) { return v.z; }" << endl
			<< "float sinTheta2(vec3 v) { return 1.0-v.z*v.z; }" << endl
			<< "float sinTheta(vec3 v) { float st2 = sinTheta2(v); if (st2 <= 0) return 0.0; else return sqrt(sinTheta2(v)); }" << endl
			<< "float tanTheta(vec3 v) { return sinTheta(v)/cosTheta(v); }" << endl
			<< endl;

		std::string vplEvalName, bsdfEvalName, lumEvalName;
		m_targetConfig.generateCode(oss, vplEvalName, bsdfEvalName, lumEvalName);

		oss << "void main() {" << endl
			<< "   /* Set up an ONB */" << endl
			<< "   vec3 N = normalize(normal);" << endl
			<< "   vec3 S = normalize(tangent - dot(tangent, N)*N);" << endl
			<< "   vec3 T = cross(N, S);" << endl
			<< endl
			<< "   /* Compute shadows */" << endl
			<< "   float d = length(lightVec);" << endl
			<< "   vec3 nLightVec = lightVec/d, absLightVec = abs(lightVec);" << endl
			<< "   float depth = max(max(absLightVec.x, absLightVec.y), absLightVec.z);" << endl
			<< "   depth = (depth-nearClip) * invClipRange - 0.001;" << endl
			<< "   float shadow = textureCube(shadowMap, nLightVec).r > depth ? 1.0 : 0.0;" << endl
			<< endl
			<< "   /* Shading */" << endl
			<< "   vec3 nCamVec = normalize(camVec);" << endl
			<< "   vec3 wo = vec3(dot(S, nLightVec)," << endl
			<< "                  dot(T, nLightVec)," << endl
			<< "                  dot(N, nLightVec));" << endl
			<< "   vec3 wi = vec3(dot(S, nCamVec)," << endl
			<< "                  dot(T, nCamVec)," << endl
			<< "                  dot(N, nCamVec));" << endl
			<< "   if (wi.z < 0)" << endl
			<< "      discard;" << endl
			<< "   vec3 vplWo = -vec3(dot(vplS, nLightVec)," << endl
			<< "                      dot(vplT, nLightVec)," << endl
			<< "                      dot(vplN, nLightVec));" << endl
			<< "   vec3 vplLo = vplPower;" << endl
			<< "   if (allowNonDiffuseVPLs)" << endl 
			<< "      vplLo *= " << vplEvalName;
			if (vpl.type == ESurfaceVPL)
				oss << "(vplUV, vplWi, vplWo);" << endl;
			else
				oss << "_dir(vplWo);" << endl;
		oss << "   if (d < minDist) d = minDist;" << endl
			<< "   gl_FragColor.rgb = vplLo * " << bsdfEvalName << "(uv, wi, wo)" << endl;
			if (vpl.type == ESurfaceVPL || (vpl.type == ELuminaireVPL 
					&& (vpl.luminaire->getType() & Luminaire::EOnSurface)))
				oss << "                      * (shadow * abs(cosTheta(wo) * cosTheta(vplWo)) / (d*d))";
			else 
				oss << "                      * (shadow * abs(cosTheta(wo)) / (d*d))";
		if (luminaire != NULL) {
			oss << endl;
			oss << "                      + " << lumEvalName << "_area(uv)"
				<< " * " << lumEvalName << "_dir(wi);" << endl;
		} else {
			oss << ";" << endl;
		}
		oss << "   gl_FragColor.a = 1.0;" << endl
			<< "}" << endl;

		program->setSource(GPUProgram::EFragmentProgram, oss.str());
		program->init();

		m_targetConfig.resolve(program);
		m_targetConfig.param_shadowMap = program->getParameterID("shadowMap", false);
		m_targetConfig.param_vplPos = program->getParameterID("vplPos", false);
		m_targetConfig.param_camPos = program->getParameterID("camPos", false);
		m_targetConfig.param_vplPower = program->getParameterID("vplPower", false);
		m_targetConfig.param_vplN = program->getParameterID("vplN", false);
		m_targetConfig.param_vplS = program->getParameterID("vplS", false);
		m_targetConfig.param_vplT = program->getParameterID("vplT", false);
		m_targetConfig.param_vplWi = program->getParameterID("vplWi", false);
		m_targetConfig.param_vplUV = program->getParameterID("vplUV", false);
		m_targetConfig.param_nearClip = program->getParameterID("nearClip", false);
		m_targetConfig.param_invClipRange = program->getParameterID("invClipRange", false);
		m_targetConfig.param_minDist = program->getParameterID("minDist", false);
		m_targetConfig.param_allowNonDiffuseVPLs = program->getParameterID("allowNonDiffuseVPLs", false);
		m_current.program = program;
		m_current.config = m_targetConfig;
		m_programs[configName] = m_current;
		program->incRef();
	}

	program->bind();
	m_shadowMap->bind(0);

	const VPLProgramConfiguration &config = m_current.config;

	program->setParameter(config.param_shadowMap, m_shadowMap);
	program->setParameter(config.param_vplPos, vpl.its.p);
	program->setParameter(config.param_camPos, camPos);
	program->setParameter(config.param_vplN, vpl.its.shFrame.n);
	program->setParameter(config.param_vplS, vpl.its.shFrame.s);
	program->setParameter(config.param_vplT, vpl.its.shFrame.t);
	if (vpl.type == ESurfaceVPL) {
		program->setParameter(config.param_vplWi, vpl.its.wi);
		program->setParameter(config.param_vplUV, vpl.its.uv);
		program->setParameter(config.param_allowNonDiffuseVPLs, m_allowNonDiffuseVPLs);
	}
	if (!m_allowNonDiffuseVPLs && vpl.type == ESurfaceVPL)
		program->setParameter(config.param_vplPower, vpl.P 
			* vpl.its.shape->getBSDF()->getDiffuseReflectance(vpl.its) * INV_PI);
	else
		program->setParameter(config.param_vplPower, vpl.P);
	program->setParameter(config.param_nearClip, m_nearClip);
	program->setParameter(config.param_invClipRange, m_invClipRange);
	program->setParameter(config.param_minDist, m_minDist);

	int textureUnitOffset = 1;
	m_targetConfig.bind(program, config, textureUnitOffset);
}

void VPLShaderManager::drawBackground(const Transform &clipToWorld, const Point &camPos) {
	if (m_backgroundProgram == NULL)
		return;
	int textureUnitOffset = 0;	
	m_backgroundProgram->bind();
	m_backgroundDependencies.recursiveBind(m_backgroundProgram, 
		m_backgroundDependencies, textureUnitOffset);
	m_backgroundProgram->setParameter("clipToWorld", clipToWorld, false);
	m_backgroundProgram->setParameter("camPos", camPos, false);
	m_renderer->blitQuad(false);
	m_backgroundProgram->unbind();
	m_backgroundDependencies.recursiveUnbind();
}

void VPLShaderManager::unbind() {
	if (m_current.program && m_current.program->isBound()) {
		m_targetConfig.unbind();
		m_current.program->unbind();
		m_shadowMap->unbind();
	}
}

void VPLShaderManager::cleanup() {
	for (std::map<std::string, ProgramAndConfiguration>::iterator it = m_programs.begin();
		it != m_programs.end(); ++it) {
		(*it).second.program->cleanup();
		(*it).second.program->decRef();
	}

	if (m_shadowMap)
		m_shadowMap->cleanup();

	if (m_backgroundProgram) {
		m_backgroundProgram->cleanup();
		m_backgroundProgram = NULL;
	}

	const std::vector<TriMesh *> meshes = m_scene->getMeshes();
	const std::vector<Luminaire *> luminaires = m_scene->getLuminaires();

	for (size_t i=0; i<meshes.size(); ++i) {
		m_renderer->unregisterGeometry(meshes[i]);
		m_renderer->unregisterShaderForResource(meshes[i]->getBSDF());
	}
	for (size_t i=0; i<luminaires.size(); ++i)
		m_renderer->unregisterShaderForResource(luminaires[i]);
	m_initialized = false;
}

MTS_IMPLEMENT_CLASS(VPLShaderManager, false, Object)
MTS_NAMESPACE_END
