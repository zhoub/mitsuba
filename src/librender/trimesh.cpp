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
#include <mitsuba/core/random.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/zstream.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/luminaire.h>

#define MTS_FILEFORMAT_HEADER 0x041C
#define MTS_FILEFORMAT_VERSION_V3 0x03

MTS_NAMESPACE_BEGIN

TriMesh::TriMesh(const std::string &name, size_t triangleCount, 
		size_t vertexCount, bool hasNormals, bool hasTexcoords, 
		bool hasVertexColors, bool flipNormals, bool faceNormals)
 	: Shape(Properties()), m_triangleCount(triangleCount),
	  m_vertexCount(vertexCount), m_flipNormals(flipNormals),
	  m_faceNormals(faceNormals) {
	m_name = name;

	m_triangles = new Triangle[m_triangleCount];
	m_positions = new Point[m_vertexCount];
	m_normals = hasNormals ? new Normal[m_vertexCount] : NULL;
	m_texcoords = hasTexcoords ? new Point2[m_vertexCount] : NULL;
	m_colors = hasVertexColors ? new Spectrum[m_vertexCount] : NULL;
	m_tangents = NULL;
}

TriMesh::TriMesh(const Properties &props) 
 : Shape(props), m_triangles(NULL), m_positions(NULL),
	m_normals(NULL), m_texcoords(NULL), m_tangents(NULL),
	m_colors(NULL) {

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

	m_triangles = NULL;
}

/* Flags used to identify available data during serialization */
enum ETriMeshFlags {
	EHasNormals      = 0x0001,
	EHasTexcoords    = 0x0002,
	EHasTangents     = 0x0004,
	EHasColors       = 0x0008,
	EFaceNormals     = 0x0010,
	ESinglePrecision = 0x1000,
	EDoublePrecision = 0x2000
};

TriMesh::TriMesh(Stream *stream, InstanceManager *manager) 
	: Shape(stream, manager), m_tangents(NULL) {
	m_name = stream->readString();
	m_aabb = AABB(stream);

	uint32_t flags = stream->readUInt();
	m_vertexCount = (size_t) stream->readULong();
	m_triangleCount = (size_t) stream->readULong();

	m_positions = new Point[m_vertexCount];
	stream->readFloatArray(reinterpret_cast<Float *>(m_positions), 
		m_vertexCount * sizeof(Point)/sizeof(Float));

	m_faceNormals = flags & EFaceNormals;

	if (flags & EHasNormals) {
		m_normals = new Normal[m_vertexCount];
		stream->readFloatArray(reinterpret_cast<Float *>(m_normals), 
			m_vertexCount * sizeof(Normal)/sizeof(Float));
	} else {
		m_normals = NULL;
	}

	if (flags & EHasTexcoords) {
		m_texcoords = new Point2[m_vertexCount];
		stream->readFloatArray(reinterpret_cast<Float *>(m_texcoords), 
			m_vertexCount * sizeof(Point2)/sizeof(Float));
	} else {
		m_texcoords = NULL;
	}

	if (flags & EHasColors) {
		m_colors = new Spectrum[m_vertexCount];
		stream->readFloatArray(reinterpret_cast<Float *>(m_colors), 
			m_vertexCount * sizeof(Spectrum)/sizeof(Float));
	} else {
		m_colors = NULL;
	}

	m_triangles = new Triangle[m_triangleCount];
	stream->readUIntArray(reinterpret_cast<uint32_t *>(m_triangles), 
		m_triangleCount * sizeof(Triangle)/sizeof(uint32_t));
	m_flipNormals = false;
	configure();
}

static void readHelper(Stream *stream, bool fileDoublePrecision, 
		Float *target, size_t count, size_t nelems) {
#if defined(SINGLE_PRECISION)
	bool hostDoublePrecision = false;
#else
	bool hostDoublePrecision = true;
#endif
	size_t size = count * nelems;
	if (fileDoublePrecision == hostDoublePrecision) {
		/* Precision matches - load directly into memory */
		stream->readFloatArray(target, size);
	} else if (fileDoublePrecision) {
		/* Double -> Single conversion */
		double *temp = new double[size];
		stream->readDoubleArray(temp, size);
		for (size_t i=0; i<size; ++i)
			target[i] = (Float) temp[i];
		delete[] temp;
	} else {
		/* Single -> Double conversion */
		float *temp = new float[size];
		stream->readSingleArray(temp, size);
		for (size_t i=0; i<size; ++i)
			target[i] = (Float) temp[i];
		delete[] temp;
	}
}

TriMesh::TriMesh(Stream *_stream) : Shape(Properties()), m_tangents(NULL) {
	ref<Stream> stream = _stream;

	if (stream->getByteOrder() != Stream::ELittleEndian) 
		Log(EError, "Tried to unserialize a shape from a stream, "
		"which was not previously set to little endian byte order!");

	if (stream->readShort() != MTS_FILEFORMAT_HEADER)
		Log(EError, "Encountered an invalid file format!");

	short version = stream->readShort();
	if (version != MTS_FILEFORMAT_VERSION_V3)
		Log(EError, "Encountered an incompatible file version!");
	stream = new ZStream(stream);

	uint32_t flags = stream->readUInt();
	m_vertexCount = (size_t) stream->readULong();
	m_triangleCount = (size_t) stream->readULong();
	
	bool fileDoublePrecision = flags & EDoublePrecision;
	m_faceNormals = flags & EFaceNormals;

	m_positions = new Point[m_vertexCount];
	readHelper(stream, fileDoublePrecision,
			reinterpret_cast<Float *>(m_positions),
			m_vertexCount, sizeof(Point)/sizeof(Float));

	if (flags & EHasNormals) {
		m_normals = new Normal[m_vertexCount];
		readHelper(stream, fileDoublePrecision, 
				reinterpret_cast<Float *>(m_normals),
				m_vertexCount, sizeof(Normal)/sizeof(Float));
	} else {
		m_normals = NULL;
	}

	if (flags & EHasTexcoords) {
		m_texcoords = new Point2[m_vertexCount];
		readHelper(stream, fileDoublePrecision,
				reinterpret_cast<Float *>(m_texcoords),
				m_vertexCount, sizeof(Point2)/sizeof(Float));
	} else {
		m_texcoords = NULL;
	}

	if (flags & EHasColors) {
		m_colors = new Spectrum[m_vertexCount];
		readHelper(stream, fileDoublePrecision, 
				reinterpret_cast<Float *>(m_colors),
				m_vertexCount, sizeof(Spectrum)/sizeof(Float));
	} else {
		m_colors = NULL;
	}

	m_triangles = new Triangle[m_triangleCount];
	stream->readUIntArray(reinterpret_cast<uint32_t *>(m_triangles), 
		m_triangleCount * sizeof(Triangle)/sizeof(uint32_t));

	m_flipNormals = false;
}

TriMesh::~TriMesh() {
	if (m_positions)
		delete[] m_positions;
	if (m_normals)
		delete[] m_normals;
	if (m_texcoords)
		delete[] m_texcoords;
	if (m_tangents)
		delete[] m_tangents;
	if (m_colors)
		delete[] m_colors;
	if (m_triangles)
		delete[] m_triangles;
}
	
std::string TriMesh::getName() const {
	return m_name;
}

AABB TriMesh::getAABB() const {
	return m_aabb;
}

Float TriMesh::pdfArea(const ShapeSamplingRecord &sRec) const {
	return m_invSurfaceArea;
}

void TriMesh::configure() {
	Shape::configure();

	if (m_areaPDF.isReady())
		return;
	
	m_aabb.reset();

	if (m_triangleCount == 0)
		Log(EError, "Encountered an empty triangle mesh!");

	/* Generate a PDF for sampling wrt. area */
	for (size_t i=0; i<m_triangleCount; i++) 
		m_areaPDF.put(m_triangles[i].surfaceArea(m_positions));
	m_surfaceArea = m_areaPDF.build();
	m_invSurfaceArea = 1.0f / m_surfaceArea;

	/* Determine the object bounds */
	for (size_t i=0; i<m_vertexCount; i++) 
		m_aabb.expandBy(m_positions[i]);
	computeNormals();

	if ((m_bsdf->getType() & BSDF::EAnisotropicMaterial
		|| m_bsdf->usesRayDifferentials()) && !m_tangents)
		computeTangentSpaceBasis();
}

Float TriMesh::getSurfaceArea() const {
	return m_surfaceArea;
}

Float TriMesh::sampleArea(ShapeSamplingRecord &sRec, const Point2 &sample) const {
	Point2 newSeed = sample;
	int index = m_areaPDF.sampleReuse(newSeed.y);
	sRec.p = m_triangles[index].sample(m_positions, m_normals, sRec.n, newSeed);
	return m_invSurfaceArea;
}

void TriMesh::computeNormals() {
	int zeroArea = 0, invalidNormals = 0;
	if (m_faceNormals) {
		if (m_normals) {
			delete[] m_normals;
			m_normals = NULL;
		}

		if (m_flipNormals) {
			/* Change the winding order */
			for (size_t i=0; i<m_triangleCount; ++i) {
				Triangle &t = m_triangles[i];
				std::swap(t.idx[0], t.idx[1]);
			}
		}
	} else {
		if (m_normals) {
			if (m_flipNormals) {
				for (size_t i=0; i<m_vertexCount; i++) 
					m_normals[i] *= -1;
			} else {
				/* Do nothing */
			}
		} else {
			m_normals = new Normal[m_vertexCount];
			memset(m_normals, 0, sizeof(Normal)*m_vertexCount);

			for (size_t i=0; i<m_triangleCount; i++) {
				const Triangle &tri = m_triangles[i];
				const Point &v0 = m_positions[tri.idx[0]];
				const Point &v1 = m_positions[tri.idx[1]];
				const Point &v2 = m_positions[tri.idx[2]];
				Normal n = Normal(cross(v1 - v0, v2 - v0));
				Float length = n.length();
				if (length != 0) {
					n /= length;
					m_normals[tri.idx[0]] += n;
					m_normals[tri.idx[1]] += n;
					m_normals[tri.idx[2]] += n;
				} else {
					zeroArea++;
				}
			}

			for (size_t i=0; i<m_vertexCount; i++) {
				Normal &n = m_normals[i];
				Float length = n.length();
				if (m_flipNormals)
					length *= -1;
				if (length != 0) {
					n /= length;
				} else {
					/* Choose some bogus value */
					invalidNormals++;
					n = Normal(1, 0, 0);
				}
			}
		}
	}

	m_flipNormals = false;
	
	if (invalidNormals > 0)
		Log(EWarn, "\"%s\": Unable to generate %i vertex normals", 
			m_name.c_str(), invalidNormals);

	if (zeroArea > 0)
		Log(EWarn, "\"%s\": Mesh contains %i zero area triangles",
			m_name.c_str(), zeroArea);
}

bool TriMesh::computeTangentSpaceBasis() {
	int zeroArea = 0, zeroNormals = 0;
	if (!m_texcoords) {
		if (m_bsdf->getType() & BSDF::EAnisotropicMaterial)
			Log(EError, "\"%s\": computeTangentSpace(): texture coordinates are required "
					"to generate tangent vectors. If you want to render with an anisotropic "
					"material, make sure that all assigned objects have texture coordinates.");
		return false;
	}

	if (m_tangents)
		Log(EError, "Tangent space vectors have already been generated!");

	m_tangents = new TangentSpace[m_vertexCount];
	memset(m_tangents, 0, sizeof(TangentSpace));

	/* No. of triangles sharing a vertex */
	uint32_t *sharers = new uint32_t[m_vertexCount];

	for (size_t i=0; i<m_vertexCount; i++) {
		m_tangents[i].dpdu = Vector(0.0f);
		m_tangents[i].dpdv = Vector(0.0f);
		if (m_normals[i].isZero()) {
			zeroNormals++;
			m_normals[i] = Normal(1.0f, 0.0f, 0.0f);
		}
		sharers[i] = 0;
	}
	
	for (size_t i=0; i<m_triangleCount; i++) {
		uint32_t idx0 = m_triangles[i].idx[0],
				 idx1 = m_triangles[i].idx[1],
				 idx2 = m_triangles[i].idx[2];
		const Point &v0 = m_positions[idx0];
		const Point &v1 = m_positions[idx1];
		const Point &v2 = m_positions[idx2];
		const Point2 &uv0 = m_texcoords[idx0];
		const Point2 &uv1 = m_texcoords[idx1];
		const Point2 &uv2 = m_texcoords[idx2];

		Vector dP1 = v1 - v0, dP2 = v2 - v0;
		Vector2 dUV1 = uv1 - uv0, dUV2 = uv2 - uv0;

		Float invDet = 1.0f, determinant = dUV1.x * dUV2.y - dUV1.y * dUV2.x;
		if (determinant != 0)
			invDet = 1.0f / determinant;

		Vector dpdu = ( dUV2.y * dP1 - dUV1.y * dP2) * invDet;
		Vector dpdv = (-dUV2.x * dP1 + dUV1.x * dP2) * invDet;

		if (dpdu.length() == 0.0f) {
			/* Recovery - required to recover from invalid geometry */
			Normal n = Normal(cross(v1 - v0, v2 - v0));
			Float length = n.length();
			if (length != 0) {
				n /= length;
				dpdu = cross(n, dpdv);
				if (dpdu.length() == 0.0f) {
					/* At least create some kind of tangent space basis 
					(fair enough for isotropic BxDFs) */
					coordinateSystem(n, dpdu, dpdv);
				}
			} else {
				zeroArea++;
			}
		}

		if (dpdv.length() == 0.0f) {
			Normal n = Normal(cross(v1 - v0, v2 - v0));
			Float length = n.length();
			if (length != 0) {
				n /= length;
				dpdv = cross(dpdu, n);
				if (dpdv.length() == 0.0f) {
					/* At least create some kind of tangent space basis 
						(fair enough for isotropic BxDFs) */
					coordinateSystem(n, dpdu, dpdv);
				}
			} else {
				zeroArea++;
			}
		}

		m_tangents[idx0].dpdu += dpdu;
		m_tangents[idx1].dpdu += dpdu;
		m_tangents[idx2].dpdu += dpdu;
		m_tangents[idx0].dpdv += dpdv;
		m_tangents[idx1].dpdv += dpdv;
		m_tangents[idx2].dpdv += dpdv;
		sharers[idx0]++; sharers[idx1]++; sharers[idx2]++;
	}

	/* Orthogonalization + Normalization pass */
	for (size_t i=0; i<m_vertexCount; i++) {
		Vector &dpdu = m_tangents[i].dpdu;
		Vector &dpdv = m_tangents[i].dpdv;

		if (dpdu.lengthSquared() == 0.0f || dpdv.lengthSquared() == 0.0f) {
			/* At least create some kind of tangent space basis 
				(fair enough for isotropic BxDFs) */
			coordinateSystem(m_normals[i], dpdu, dpdv);
		} else {
			if (sharers[i] > 0) {
				dpdu /= (Float) sharers[i];
				dpdv /= (Float) sharers[i];
			}
		}
	}
	delete[] sharers;

	if (zeroArea > 0 || zeroNormals > 0)
		Log(EWarn, "\"%s\": computeTangentSpace(): Mesh contains invalid "
			"geometry: %i zero area triangles and %i zero normals found!", 
			m_name.c_str(), zeroArea, zeroNormals);
	return true;
}

ref<TriMesh> TriMesh::createTriMesh() {
	return this;
}

void TriMesh::serialize(Stream *stream, InstanceManager *manager) const {
	Shape::serialize(stream, manager);
	uint32_t flags = 0;
	if (m_normals)
		flags |= EHasNormals;
	if (m_texcoords)
		flags |= EHasTexcoords;
	if (m_colors)
		flags |= EHasColors;
	if (m_faceNormals)
		flags |= EFaceNormals;
	stream->writeString(m_name);
	m_aabb.serialize(stream);
	stream->writeUInt(flags);
	stream->writeULong(m_vertexCount);
	stream->writeULong(m_triangleCount);

	stream->writeFloatArray(reinterpret_cast<Float *>(m_positions), 
		m_vertexCount * sizeof(Point)/sizeof(Float));
	if (m_normals)
		stream->writeFloatArray(reinterpret_cast<Float *>(m_normals), 
			m_vertexCount * sizeof(Normal)/sizeof(Float));
	if (m_texcoords)
		stream->writeFloatArray(reinterpret_cast<Float *>(m_texcoords), 
			m_vertexCount * sizeof(Point2)/sizeof(Float));
	if (m_colors)
		stream->writeFloatArray(reinterpret_cast<Float *>(m_colors), 
			m_vertexCount * sizeof(Spectrum)/sizeof(Float));
	stream->writeUIntArray(reinterpret_cast<uint32_t *>(m_triangles), 
		m_triangleCount * sizeof(Triangle)/sizeof(uint32_t));
}

void TriMesh::writeOBJ(const fs::path &path) const {
	fs::ofstream os(path);
	os << "o " << m_name << endl;
	for (size_t i=0; i<m_vertexCount; ++i) {
		os << "v " 
			<< m_positions[i].x << " "
			<< m_positions[i].y << " "
			<< m_positions[i].z << endl;
	}

	if (m_normals) {
		for (size_t i=0; i<m_vertexCount; ++i) {
			os << "vn " 
				<< m_normals[i].x << " "
				<< m_normals[i].y << " "
				<< m_normals[i].z << endl;
		}
	}

	if (m_normals) {
		for (size_t i=0; i<m_triangleCount; ++i) {
			os << "f " 
				<< m_triangles[i].idx[0] + 1 << "//" 
				<< m_triangles[i].idx[0] + 1 << " "
				<< m_triangles[i].idx[1] + 1 << "//" 
				<< m_triangles[i].idx[1] + 1 << " "
				<< m_triangles[i].idx[2] + 1 << "//" 
				<< m_triangles[i].idx[2] + 1 << endl;
		}
	} else {
		for (size_t i=0; i<m_triangleCount; ++i) {
			os << "f " 
				<< m_triangles[i].idx[0] + 1 << " "
				<< m_triangles[i].idx[1] + 1 << " "
				<< m_triangles[i].idx[2] + 1 << endl;
		}
	}

	os.close();
}

void TriMesh::serialize(Stream *_stream) const {
	ref<Stream> stream = _stream;

	if (stream->getByteOrder() != Stream::ELittleEndian) 
		Log(EError, "Tried to unserialize a shape from a stream, "
			"which was not previously set to little endian byte order!");

	stream->writeShort(MTS_FILEFORMAT_HEADER);
	stream->writeShort(MTS_FILEFORMAT_VERSION_V3);
	stream = new ZStream(stream);

#if defined(SINGLE_PRECISION)
	uint32_t flags = ESinglePrecision;
#else
	uint32_t flags = EDoublePrecision;
#endif

	if (m_normals)
		flags |= EHasNormals;
	if (m_texcoords)
		flags |= EHasTexcoords;
	if (m_colors)
		flags |= EHasColors;
	if (m_faceNormals)
		flags |= EFaceNormals;

	stream->writeUInt(flags);
	stream->writeULong(m_vertexCount);
	stream->writeULong(m_triangleCount);

	stream->writeFloatArray(reinterpret_cast<Float *>(m_positions), 
		m_vertexCount * sizeof(Point)/sizeof(Float));
	if (m_normals)
		stream->writeFloatArray(reinterpret_cast<Float *>(m_normals), 
			m_vertexCount * sizeof(Normal)/sizeof(Float));
	if (m_texcoords)
		stream->writeFloatArray(reinterpret_cast<Float *>(m_texcoords), 
			m_vertexCount * sizeof(Point2)/sizeof(Float));
	if (m_colors)
		stream->writeFloatArray(reinterpret_cast<Float *>(m_colors), 
			m_vertexCount * sizeof(Spectrum)/sizeof(Float));
	stream->writeUIntArray(reinterpret_cast<uint32_t *>(m_triangles), 
		m_triangleCount * sizeof(Triangle)/sizeof(uint32_t));
}

std::string TriMesh::toString() const {
	std::ostringstream oss;
	oss << getClass()->getName() << "[" << endl
		<< "  name = \"" << m_name<< "\"," << endl
		<< "  triangleCount = " << m_triangleCount << "," << endl
		<< "  vertexCount = " << m_vertexCount << "," << endl
		<< "  faceNormals = " << (m_faceNormals ? "true" : "false") << "," << endl
		<< "  hasNormals = " << (m_normals ? "true" : "false") << "," << endl
		<< "  hasTexcoords = " << (m_texcoords ? "true" : "false") << "," << endl
		<< "  hasTangents = " << (m_tangents ? "true" : "false") << "," << endl
		<< "  hasColors = " << (m_colors ? "true" : "false") << "," << endl
		<< "  surfaceArea = " << m_surfaceArea << "," << endl
		<< "  aabb = " << m_aabb.toString() << "," << endl
		<< "  bsdf = " << indent(m_bsdf.toString()) << "," << endl
		<< "  subsurface = " << indent(m_subsurface.toString()) << "," << endl
		<< "  luminaire = " << indent(m_luminaire.toString()) << endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS_S(TriMesh, false, Shape)
MTS_NAMESPACE_END
