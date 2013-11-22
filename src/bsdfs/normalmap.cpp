/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

class NormalMap : public BSDF {
public:
	NormalMap(const Properties &props) : BSDF(props) { }

	NormalMap(Stream *stream, InstanceManager *manager)
			: BSDF(stream, manager) {
		m_nested = static_cast<BSDF *>(manager->getInstance(stream));
		m_normals = static_cast<Texture *>(manager->getInstance(stream));
		configure();
	}

	void configure() {
		if (!m_nested)
			Log(EError, "A child BSDF instance is required");
		if (!m_normals)
			Log(EError, "A normal map texture must be specified");

		m_components.clear();
		for (int i=0; i<m_nested->getComponentCount(); ++i)
			m_components.push_back(m_nested->getType(i) | ESpatiallyVarying | EAnisotropic);

		m_usesRayDifferentials = true;

		BSDF::configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		manager->serialize(stream, m_nested.get());
		manager->serialize(stream, m_normals.get());
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
			if (m_nested != NULL)
				Log(EError, "Only a single nested BSDF can be added!");
			m_nested = static_cast<BSDF *>(child);
		} else if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (m_normals != NULL)
				Log(EError, "Only a single normals texture can be specified!");
			m_normals = static_cast<Texture *>(child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	Frame getFrame(const Intersection &its) const {
		Normal n;
		m_normals->eval(its, false).toLinearRGB(n.x, n.y, n.z);
		for (int i=0; i<3; ++i)
			n[i] = 2 * n[i] - 1;

		Frame result;
		result.n = normalize(its.shFrame.toWorld(n));
		result.s = normalize(its.dpdu - result.n
			* dot(result.n, its.dpdu));
		result.t = cross(result.n, result.s);

		if (dot(result.n, its.geoFrame.n) < 0)
			result.n *= -1;

		return result;
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		const Intersection& its = bRec.its;
		Intersection perturbed(its);
		perturbed.shFrame = getFrame(its);

		BSDFSamplingRecord perturbedQuery(perturbed,
			perturbed.toLocal(its.toWorld(bRec.wi)),
			perturbed.toLocal(its.toWorld(bRec.wo)), bRec.mode);

		if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(perturbedQuery.wo) <= 0)
			return Spectrum(0.0f);

		perturbedQuery.sampler = bRec.sampler;
		perturbedQuery.typeMask = bRec.typeMask;
		perturbedQuery.component = bRec.component;

		return m_nested->eval(perturbedQuery, measure);
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		const Intersection& its = bRec.its;
		Intersection perturbed(its);
		perturbed.shFrame = getFrame(its);

		BSDFSamplingRecord perturbedQuery(perturbed,
			perturbed.toLocal(its.toWorld(bRec.wi)),
			perturbed.toLocal(its.toWorld(bRec.wo)), bRec.mode);
		if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(perturbedQuery.wo) <= 0)
			return 0;
		perturbedQuery.mode = bRec.mode;
		perturbedQuery.sampler = bRec.sampler;
		perturbedQuery.typeMask = bRec.typeMask;
		perturbedQuery.component = bRec.component;
		return m_nested->pdf(perturbedQuery, measure);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		const Intersection& its = bRec.its;
		Intersection perturbed(its);
		perturbed.shFrame = getFrame(its);

		BSDFSamplingRecord perturbedQuery(perturbed, bRec.sampler, bRec.mode);
		perturbedQuery.wi = perturbed.toLocal(its.toWorld(bRec.wi));
		perturbedQuery.sampler = bRec.sampler;
		perturbedQuery.typeMask = bRec.typeMask;
		perturbedQuery.component = bRec.component;
		Spectrum result = m_nested->sample(perturbedQuery, sample);
		if (!result.isZero()) {
			bRec.sampledComponent = perturbedQuery.sampledComponent;
			bRec.sampledType = perturbedQuery.sampledType;
			bRec.wo = its.toLocal(perturbed.toWorld(perturbedQuery.wo));
			bRec.eta = perturbedQuery.eta;
			if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(perturbedQuery.wo) <= 0)
				return Spectrum(0.0f);
		}
		return result;
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		const Intersection& its = bRec.its;
		Intersection perturbed(its);
		perturbed.shFrame = getFrame(its);

		BSDFSamplingRecord perturbedQuery(perturbed, bRec.sampler, bRec.mode);
		perturbedQuery.wi = perturbed.toLocal(its.toWorld(bRec.wi));
		perturbedQuery.typeMask = bRec.typeMask;
		perturbedQuery.component = bRec.component;
		Spectrum result = m_nested->sample(perturbedQuery, pdf, sample);

		if (!result.isZero()) {
			bRec.sampledComponent = perturbedQuery.sampledComponent;
			bRec.sampledType = perturbedQuery.sampledType;
			bRec.wo = its.toLocal(perturbed.toWorld(perturbedQuery.wo));
			bRec.eta = perturbedQuery.eta;
			if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(perturbedQuery.wo) <= 0)
				return Spectrum(0.0f);
		}

		return result;
	}

	Float getRoughness(const Intersection &its, int component) const {
		return m_nested->getRoughness(its, component);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "NormalMap[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  normals = " << indent(m_normals->toString()) << endl
			<< "  nested = " << indent(m_nested->toString()) << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
protected:
	ref<Texture> m_normals;
	ref<BSDF> m_nested;
};

// ================ Hardware shader implementation ================

/**
 * This is a quite approximate version of the bump map model -- it likely
 * won't match the reference exactly, but it should be good enough for
 * preview purposes
 */
class NormalMapShader : public Shader {
public:
	NormalMapShader(Renderer *renderer, const BSDF *nested, const Texture *normals)
		: Shader(renderer, EBSDFShader), m_nested(nested), m_normals(normals) {
		m_nestedShader = renderer->registerShaderForResource(m_nested.get());
		m_normalShader = renderer->registerShaderForResource(m_normals.get());
	}

	bool isComplete() const {
		return m_nestedShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_nested.get());
		renderer->unregisterShaderForResource(m_normals.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_nestedShader.get());
		deps.push_back(m_normalShader.get());
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    vec3 n = normalize(2.0*" << depNames[1] << "(uv) - vec3(1.0));" << endl
			<< "    vec3 s = normalize(vec3(1.0-n.x*n.x, -n.x*n.y, -n.x*n.z)); " << endl
			<< "    vec3 t = cross(s, n);" << endl
			<< "    wi = vec3(dot(wi, s), dot(wi, t), dot(wi, n));" << endl
			<< "    wo = vec3(dot(wo, s), dot(wo, t), dot(wo, n));" << endl
			<< "    return " << depNames[0] << "(uv, wi, wo);" << endl
			<< "}" << endl
			<< endl
		    << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    vec3 n = normalize(2.0*" << depNames[1] << "(uv) - vec3(1.0));" << endl
			<< "    vec3 s = normalize(vec3(1.0-n.x*n.x, -n.x*n.y, -n.x*n.z)); " << endl
			<< "    vec3 t = cross(s, n);" << endl
			<< "    wi = vec3(dot(wi, s), dot(wi, t), dot(wi, n));" << endl
			<< "    wo = vec3(dot(wo, s), dot(wo, t), dot(wo, n));" << endl
			<< "    return " << depNames[0] << "_diffuse(uv, wi, wo);" << endl
			<< "}" << endl
			<< endl;
	}

	MTS_DECLARE_CLASS()
private:
	ref<const BSDF> m_nested;
	ref<const Texture> m_normals;
	ref<Shader> m_nestedShader;
	ref<Shader> m_normalShader;
};

Shader *NormalMap::createShader(Renderer *renderer) const {
	return new NormalMapShader(renderer, m_nested.get(), m_normals.get());
}

MTS_IMPLEMENT_CLASS(NormalMapShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(NormalMap, false, BSDF)
MTS_EXPORT_PLUGIN(NormalMap, "Smooth dielectric coating");
MTS_NAMESPACE_END
