/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/luminaire.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <set>

MTS_NAMESPACE_BEGIN

/**
 * Wavefront OBJ triangle mesh loader
 */
class WavefrontOBJ : public Shape {
public:
	struct OBJTriangle {
		unsigned int p[3];
		unsigned int n[3];
		unsigned int uv[3];
	};

	WavefrontOBJ(const Properties &props) : Shape(props) {
		FileResolver *fResolver = Thread::getThread()->getFileResolver();
		fs::path path = fResolver->resolve(props.getString("filename"));
		m_name = path.stem();
	
		/* By default, any existing normals will be used for
		   rendering. If no normals are found, Mitsuba will
		   automatically generate smooth vertex normals. 
		   Setting the 'faceNormals' parameter instead forces
		   the use of face normals, which will result in a faceted
		   appearance.
		*/
		m_faceNormals = props.getBoolean("faceNormals", false);

		/* Causes all normals to be flipped */
		m_flipNormals = props.getBoolean("flipNormals", false);

		/* Load the geometry */
		Log(EInfo, "Loading geometry from \"%s\" ..", path.leaf().c_str());
		fs::ifstream is(path);
		if (is.bad() || is.fail())
			Log(EError, "Geometry file '%s' not found!", path.file_string().c_str());

		std::string buf;
		std::vector<Point> vertices;
		std::vector<Normal> normals;
		std::vector<Point2> texcoords;
		std::vector<OBJTriangle> triangles;
		bool hasNormals = false, hasTexcoords = false;
		bool firstVertex = true;
		BSDF *currentMaterial = NULL;
		std::string name = m_name;
		std::set<std::string> geomNames;
		int geomIdx = 0;

		while (is >> buf) {
			if (buf == "v") {
				/* Parse + transform vertices */
				Point p;
				is >> p.x >> p.y >> p.z;
				vertices.push_back(m_objectToWorld(p));
				if (firstVertex) {
					if (triangles.size() > 0) {
						if (geomNames.find(name) != geomNames.end())
							/// make sure that we have unique names
							name = formatString("%s_%i", m_name.c_str(), geomIdx);
						generateGeometry(name, vertices, normals, texcoords, 
							triangles, hasNormals, hasTexcoords, currentMaterial);
						triangles.clear();
						geomNames.insert(name);
						geomIdx++;
					}
					hasNormals = false;
					hasTexcoords = false;
					firstVertex = false;
				}
			} else if (buf == "vn") {
				Normal n;
				is >> n.x >> n.y >> n.z;
				if (!n.isZero())
					normals.push_back(normalize(m_objectToWorld(n)));
				else
					normals.push_back(n);
				hasNormals = true;
			} else if (buf == "g") {
				std::string line;
				std::getline(is, line);
				if (line.length() > 2) { 
					name = trim(line.substr(1, line.length()-1));
					Log(EInfo, "Loading geometry \"%s\"", name.c_str());
				}
			} else if (buf == "usemtl") {
				std::string line;
				std::getline(is, line);
				std::string materialName = trim(line.substr(1, line.length()-1));
				if (m_materials.find(materialName) != m_materials.end())
					currentMaterial = m_materials[materialName];
				else
					currentMaterial = NULL;
			} else if (buf == "mtllib") {
				std::string line;
				std::getline(is, line);
				ref<FileResolver> frClone = fResolver->clone();
				frClone->addPath(fs::complete(path).parent_path());
				fs::path mtlName = frClone->resolve(trim(line.substr(1, line.length()-1)));
				if (fs::exists(mtlName))
					parseMaterials(mtlName);
				else
					Log(EWarn, "Could not find referenced material library '%s'", 
						mtlName.file_string().c_str());
			} else if (buf == "vt") {
				std::string line;
				Float u, v, w;
				std::getline(is, line);
				std::istringstream iss(line);
				iss >> u >> v >> w;
				texcoords.push_back(Point2(u, v));
				hasTexcoords = true;
			} else if (buf == "f") {
				std::string line, tmp;
				std::getline(is, line);
				std::istringstream iss(line);
				firstVertex = true;
				OBJTriangle t;
				iss >> tmp; parse(t, 0, tmp);
				iss >> tmp; parse(t, 1, tmp);
				iss >> tmp; parse(t, 2, tmp);
				triangles.push_back(t);
				if (iss >> tmp) {
					parse(t, 1, tmp);
					std::swap(t.p[0], t.p[1]);
					std::swap(t.uv[0], t.uv[1]);
					std::swap(t.n[0], t.n[1]);
					triangles.push_back(t);
				}
			} else {
				/* Ignore */
				std::string line;
				std::getline(is, line);
			}
		}

		generateGeometry(name, vertices, normals, texcoords, 
			triangles, hasNormals, hasTexcoords, currentMaterial);
	}

	WavefrontOBJ(Stream *stream, InstanceManager *manager) : Shape(stream, manager) {
		unsigned int meshCount = stream->readUInt();
		m_meshes.resize(meshCount);
	
		for (unsigned int i=0; i<meshCount; ++i) {
			m_meshes[i] = static_cast<TriMesh *>(manager->getInstance(stream));
			m_meshes[i]->incRef();
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Shape::serialize(stream, manager);

		stream->writeUInt((unsigned int) m_meshes.size());
		for (size_t i=0; i<m_meshes.size(); ++i)
			manager->serialize(stream, m_meshes[i]);
	}

	void parse(OBJTriangle &t, int i, const std::string &str) {
		std::vector<std::string> tokens = tokenize(str, "/");
		if (tokens.size() == 1) {
			t.p[i] = atoi(tokens[0].c_str())-1;
		} else if (tokens.size() == 2) {
			if (str.find("//") == std::string::npos) {
				t.p[i]  = atoi(tokens[0].c_str())-1;
				t.uv[i] = atoi(tokens[1].c_str())-1;
			} else {
				t.p[i] = atoi(tokens[0].c_str())-1;
				t.n[i] = atoi(tokens[1].c_str())-1;
			}
		} else if (tokens.size() == 3) {
			t.p[i] = atoi(tokens[0].c_str())-1;
			t.uv[i] = atoi(tokens[1].c_str())-1;
			t.n[i] = atoi(tokens[2].c_str())-1;
		} else {
			Log(EError, "Invalid OBJ face format!");
		}
	}

	void parseMaterials(const fs::path &mtlPath) {
		Log(EInfo, "Loading OBJ materials from \"%s\" ..", mtlPath.filename().c_str());
		fs::ifstream is(mtlPath);
		if (is.bad() || is.fail())
			Log(EError, "Unexpected I/O error while accessing material file '%s'!", 
				mtlPath.file_string().c_str());
		std::string buf;
		std::string mtlName;
		Spectrum diffuse;

		while (is >> buf) {
			if (buf == "newmtl") {
				if (mtlName != "") 
					addMaterial(mtlName, diffuse);

				std::string line, tmp;
				std::getline(is, line);
				mtlName = trim(line.substr(1, line.length()-1));
			} else if (buf == "Kd") {
				Float r, g, b;
				is >> r >> g >> b;
				diffuse.fromSRGB(r, g, b);
			} else {
				/* Ignore */
				std::string line;
				std::getline(is, line);
			}
		}
		addMaterial(mtlName, diffuse);
	}

	void addMaterial(const std::string &name, const Spectrum &diffuse) {
		Properties props("lambertian");
		props.setSpectrum("reflectance", diffuse);
		props.setID(name);
		BSDF *bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(BSDF::m_theClass, props));
		bsdf->incRef();
		m_materials[name] = bsdf;
	}

	struct Vertex {
		Point p;
		Normal n;
		Point2 uv;
	};

	/// For using vertices as keys in an associative structure
	struct vertex_key_order : public 
		std::binary_function<Vertex, Vertex, bool> {
	public:
		bool operator()(const Vertex &v1, const Vertex &v2) const {
			if (v1.p.x < v2.p.x) return true;
			else if (v1.p.x > v2.p.x) return false;
			if (v1.p.y < v2.p.y) return true;
			else if (v1.p.y > v2.p.y) return false;
			if (v1.p.z < v2.p.z) return true;
			else if (v1.p.z > v2.p.z) return false;
			if (v1.n.x < v2.n.x) return true;
			else if (v1.n.x > v2.n.x) return false;
			if (v1.n.y < v2.n.y) return true;
			else if (v1.n.y > v2.n.y) return false;
			if (v1.n.z < v2.n.z) return true;
			else if (v1.n.z > v2.n.z) return false;
			if (v1.uv.x < v2.uv.x) return true;
			else if (v1.uv.x > v2.uv.x) return false;
			if (v1.uv.y < v2.uv.y) return true;
			else if (v1.uv.y > v2.uv.y) return false;
			return false;
		}
	};

	void generateGeometry(const std::string &name,
			const std::vector<Point> &vertices,
			const std::vector<Normal> &normals,
			const std::vector<Point2> &texcoords,
			const std::vector<OBJTriangle> &triangles,
			bool hasNormals, bool hasTexcoords,
			BSDF *currentMaterial) {
		if (triangles.size() == 0)
			return;

		std::map<Vertex, int, vertex_key_order> vertexMap;
		std::vector<Vertex> vertexBuffer;
		size_t numMerged = 0;

		/* Collapse the mesh into a more usable form */
		Triangle *triangleArray = new Triangle[triangles.size()];
		for (unsigned int i=0; i<triangles.size(); i++) {
			Triangle tri;
			for (unsigned int j=0; j<3; j++) {
				unsigned int vertexId = triangles[i].p[j];
				unsigned int normalId = triangles[i].n[j];
				unsigned int uvId = triangles[i].uv[j];
				int key;

				Vertex vertex;
				vertex.p = vertices.at(vertexId);
				if (hasNormals)
					vertex.n = normals.at(normalId);
				else
					vertex.n = Normal(0.0f);
				if (hasTexcoords)
					vertex.uv = texcoords.at(uvId);
				else
					vertex.uv = Point2(0.0f);

				if (vertexMap.find(vertex) != vertexMap.end()) {
					key = vertexMap[vertex];
					numMerged++;
				} else {
					key = (int) vertexBuffer.size();
					vertexMap[vertex] = (int) key;
					vertexBuffer.push_back(vertex);
				}

				tri.idx[j] = key;
			}
			triangleArray[i] = tri;
		}

		ref<TriMesh> mesh = new TriMesh(name, m_worldToObject,
			triangles.size(), vertexBuffer.size(),
			hasNormals, hasTexcoords, false,
			m_flipNormals, m_faceNormals);

		std::copy(triangleArray, triangleArray+triangles.size(), mesh->getTriangles());

		Point    *target_positions = mesh->getVertexPositions();
		Normal   *target_normals   = mesh->getVertexNormals();
		Point2   *target_texcoords = mesh->getVertexTexcoords();

		for (size_t i=0; i<vertexBuffer.size(); i++) {
			*target_positions++ = vertexBuffer[i].p;
			if (hasNormals)
				*target_normals++ = vertexBuffer[i].n;
			if (hasTexcoords)
				*target_texcoords++ = vertexBuffer[i].uv;
		}

		mesh->incRef();
		mesh->configure();
		if (currentMaterial)
			mesh->addChild("", currentMaterial);
		m_meshes.push_back(mesh);
		Log(EInfo, "%s: Loaded " SIZE_T_FMT " triangles, " SIZE_T_FMT 
			" vertices (merged " SIZE_T_FMT " vertices).", name.c_str(),
			triangles.size(), vertexBuffer.size(), numMerged);
	}

	virtual ~WavefrontOBJ() {
		for (size_t i=0; i<m_meshes.size(); ++i)
			m_meshes[i]->decRef();
		for (std::map<std::string, BSDF *>::iterator it = m_materials.begin();
			it != m_materials.end(); ++it) {
			(*it).second->decRef();
		}
	}
	
	void configure() {
		Shape::configure();

		m_aabb.reset();
		for (size_t i=0; i<m_meshes.size(); ++i) {
			m_meshes[i]->configure();
			m_aabb.expandBy(m_meshes[i]->getAABB());
		}
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		const Class *cClass = child->getClass();
		if (cClass->derivesFrom(BSDF::m_theClass)) {
			m_bsdf = static_cast<BSDF *>(child);
			for (size_t i=0; i<m_meshes.size(); ++i) 
				m_meshes[i]->addChild(name, child);
			Assert(m_meshes.size() > 0);
			m_bsdf->setParent(NULL);
		} else if (cClass->derivesFrom(Luminaire::m_theClass)) {
			Assert(m_luminaire == NULL && m_meshes.size() == 1);
			m_luminaire = static_cast<Luminaire *>(child);
			for (size_t i=0; i<m_meshes.size(); ++i) {
				child->setParent(m_meshes[i]);
				m_meshes[i]->addChild(name, child);
			}
		} else if (cClass->derivesFrom(Subsurface::m_theClass)) {
			Assert(m_subsurface == NULL);
			m_subsurface = static_cast<Subsurface *>(child);
			for (size_t i=0; i<m_meshes.size(); ++i) { 
				child->setParent(m_meshes[i]);
				m_meshes[i]->addChild(name, child);
			}
		} else {
			Shape::addChild(name, child);
		}
	}

	bool isCompound() const {
		return true;
	}

	Shape *getElement(int index) {
		if (index >= (int) m_meshes.size())
			return NULL;
		Shape *shape = m_meshes[index];
		BSDF *bsdf = shape->getBSDF();
		Luminaire *luminaire = shape->getLuminaire();
		Subsurface *subsurface = shape->getSubsurface();
		if (bsdf)
			bsdf->setParent(shape);
		if (luminaire)
			luminaire->setParent(shape);
		if (subsurface)
			subsurface->setParent(shape);
		return shape;
	}

	MTS_DECLARE_CLASS()
private:
	std::vector<TriMesh *> m_meshes;
	std::map<std::string, BSDF *> m_materials;
	bool m_flipNormals, m_faceNormals;
};

MTS_IMPLEMENT_CLASS_S(WavefrontOBJ, false, Shape)
MTS_EXPORT_PLUGIN(WavefrontOBJ, "OBJ triangle mesh loader");
MTS_NAMESPACE_END
