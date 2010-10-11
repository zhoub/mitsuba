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
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/timer.h>
#include <ply/ply_parser.hpp>
#include <tr1/functional>

using namespace std::tr1::placeholders;

MTS_NAMESPACE_BEGIN

/**
 * PLY mesh loader using libply by Ares Lagae
 * (http://people.cs.kuleuven.be/~ares.lagae/libply/)
 */
class PLYLoader : public TriMesh {
public:
	PLYLoader(const Properties &props) : TriMesh(props) {
		fs::path filePath = Thread::getThread()->getFileResolver()->resolve(
			props.getString("filename"));
		m_name = filePath.stem();

		/* Determines whether vertex colors should be 
		   treated as linear RGB or sRGB. */
		m_sRGB = props.getBoolean("srgb", true);

		/* Load the geometry */
		Log(EInfo, "Loading geometry from \"%s\" ..", filePath.leaf().c_str());
		m_triangleCount = m_vertexCount = 0;
		m_vertexCtr = m_triangleCtr = m_triangleIdxCtr = 0;
		m_normal = Normal(0.0f);
		m_uv = Point2(0.0f);
		m_hasNormals = false;
		m_hasTexCoords = false;
		memset(&m_triangle, 0, sizeof(Triangle));
		loadPLY(filePath);
		if (m_triangleCount == 0 || m_vertexCount == 0)
			Log(EError, "Unable to load \"%s\" (no triangles or vertices found)!");

		Assert(m_triangleCtr == m_triangleCount);
		Assert(m_vertexCtr == m_vertexCount);

		calculateTangentSpaceBasis(m_hasNormals, m_hasTexCoords, true);
	}

	PLYLoader(Stream *stream, InstanceManager *manager) : TriMesh(stream, manager) { }

	void loadPLY(const fs::path &path);

	void info_callback(const std::string& filename, std::size_t line_number,
			const std::string& message) {
		Log(EInfo, "\"%s\" [line %i] info: %s", filename.c_str(), line_number,
			message.c_str());
	}

	void warning_callback(const std::string& filename, std::size_t line_number,
			const std::string& message) {
		Log(EWarn, "\"%s\" [line %i] warning: %s", filename.c_str(), line_number,
			message.c_str());
	}

	void error_callback(const std::string& filename, std::size_t line_number,
			const std::string& message) {
		Log(EError, "\"%s\" [line %i] error: %s", filename.c_str(), line_number,
			message.c_str());
	}

	template<typename ValueType> std::tr1::function <void (ValueType)> 
		scalar_property_definition_callback(const std::string& element_name, 
		const std::string& property_name);

	template<typename SizeType, typename IndexType> std::tr1::tuple<std::tr1::function<void (SizeType)>, 
		std::tr1::function<void (IndexType)>, std::tr1::function<void ()> > 
		list_property_definition_callback(const std::string& element_name,
		const std::string& property_name);

	std::tr1::tuple<std::tr1::function<void()>, std::tr1::function<void()> > 
		element_definition_callback(const std::string& element_name, std::size_t count) {
		if (element_name == "vertex") {
			m_vertexCount = count;
			m_vertexBuffer = new Vertex[count];
			return std::tr1::tuple<std::tr1::function<void()>,
				std::tr1::function<void()> >(
					std::tr1::bind(&PLYLoader::vertex_begin_callback, this),
					std::tr1::bind(&PLYLoader::vertex_end_callback, this)
			);
		} else if (element_name == "face") {
			m_triangleCount = count;
			m_triangles = new Triangle[m_triangleCount];
			return std::tr1::tuple<std::tr1::function<void()>,
				std::tr1::function<void()> >(
				std::tr1::bind(&PLYLoader::face_begin_callback, this),
				std::tr1::bind(&PLYLoader::face_end_callback, this)
			);
		} else {
			return
				std::tr1::tuple<std::tr1::function<void()>,
				std::tr1::function<void()> >(0, 0);
		}
	}

	void vertex_x_callback(ply::float32 x) { m_position.x = x; }
	void vertex_y_callback(ply::float32 y) { m_position.y = y; }
	void vertex_z_callback(ply::float32 z) { m_position.z = z; }
	void normal_x_callback(ply::float32 x) { m_normal.x = x; }
	void normal_y_callback(ply::float32 y) { m_normal.y = y; }
	void normal_z_callback(ply::float32 z) { m_normal.z = z; }
	void texcoord_u_callback(ply::float32 x) { m_uv.x = x; }
	void texcoord_v_callback(ply::float32 y) { m_uv.y = y; }

	inline Float fromSRGBComponent(Float value) {
		if (value <= (Float) 0.04045)
			return value / (Float) 12.92;
		return std::pow((value + (Float) 0.055)
			/ (Float) (1.0 + 0.055), (Float) 2.4);
	}

	void vertex_begin_callback() { }
	void vertex_end_callback() {
		m_vertexBuffer[m_vertexCtr].p = m_objectToWorld(m_position);
		m_vertexBuffer[m_vertexCtr].n = m_objectToWorld(m_normal);
		m_vertexBuffer[m_vertexCtr].uv = m_uv;
		m_vertexBuffer[m_vertexCtr].dpdu = Vector(0.0f);
		m_vertexBuffer[m_vertexCtr].dpdv = Vector(0.0f);
#if defined(MTS_HAS_VERTEX_COLORS)
		if (m_sRGB) {
			m_red = fromSRGBComponent(m_red);
			m_green = fromSRGBComponent(m_green);
			m_blue = fromSRGBComponent(m_blue);
		}
		m_vertexBuffer[m_vertexCtr].color[0] = m_red;
		m_vertexBuffer[m_vertexCtr].color[1] = m_green;
		m_vertexBuffer[m_vertexCtr].color[2] = m_blue;
#endif
		m_vertexCtr++;
	}

	void red_callback_uint8(ply::uint8 r) { m_red = r / 255.0f; }
	void green_callback_uint8(ply::uint8 g) { m_green = g / 255.0f; }
	void blue_callback_uint8(ply::uint8 b) { m_blue = b / 255.0f; }

	void red_callback(ply::float32 r) { m_red = r; }
	void green_callback(ply::float32 g) { m_green = g; }
	void blue_callback(ply::float32 b) { m_blue = b; }

	void face_begin_callback() { }
	void face_end_callback() { }

	void face_vertex_indices_begin_uint8(ply::uint8 size) { 
		if (size != 3)
			Log(EError, "Only triangle PLY meshes are supported for now.");
		m_triangleIdxCtr = 0;
	}

	void face_vertex_indices_begin_uint32(ply::uint32 size) { 
		if (size != 3)
			Log(EError, "Only triangle PLY meshes are supported for now.");
		m_triangleIdxCtr = 0;
	}

	void face_vertex_indices_element_int32(ply::int32 element) {
		Assert(m_triangleIdxCtr < 3);
		Assert((size_t) element < m_vertexCount);
		m_triangle.idx[m_triangleIdxCtr++] = element;
	}

	void face_vertex_indices_element_uint32(ply::uint32 element) {
		Assert(m_triangleIdxCtr < 3);
		Assert((size_t) element < m_vertexCount);
		m_triangle.idx[m_triangleIdxCtr++] = element;
	}

	void face_vertex_indices_end() { 
		Assert(m_triangleIdxCtr == 3);
		m_triangles[m_triangleCtr++] = m_triangle;
	}

	MTS_DECLARE_CLASS()
private:
	Point m_position;
	Normal m_normal;
	Float m_red, m_green, m_blue;
	size_t m_vertexCtr, m_triangleCtr, m_triangleIdxCtr;
	Triangle m_triangle;
	bool m_hasNormals, m_hasTexCoords;
	Point2 m_uv;
	bool m_sRGB;
};
	
template<> std::tr1::function <void (ply::float32)> 
	PLYLoader::scalar_property_definition_callback(const std::string& element_name, 
	const std::string& property_name) {
	if (element_name == "vertex") {
		if (property_name == "x") {
			return std::tr1::bind(&PLYLoader::vertex_x_callback, this,  _1);
		} else if (property_name == "y") {
			return std::tr1::bind(&PLYLoader::vertex_y_callback, this,  _1);
		} else if (property_name == "z") {
			return std::tr1::bind(&PLYLoader::vertex_z_callback, this,  _1);
		} else if (property_name == "nx") {
			m_hasNormals = true;
			return std::tr1::bind(&PLYLoader::normal_x_callback, this,  _1);
		} else if (property_name == "ny") {
			return std::tr1::bind(&PLYLoader::normal_y_callback, this,  _1);
		} else if (property_name == "nz") {
			return std::tr1::bind(&PLYLoader::normal_z_callback, this,  _1);
		} else if (property_name == "u") {
			m_hasTexCoords = true;
			return std::tr1::bind(&PLYLoader::texcoord_u_callback, this,  _1);
		} else if (property_name == "v") {
			return std::tr1::bind(&PLYLoader::texcoord_v_callback, this,  _1);
		} else if (property_name == "diffuse_red" || property_name == "red") {
			return std::tr1::bind(&PLYLoader::red_callback, this,  _1);
		} else if (property_name == "diffuse_green" || property_name == "green") {
			return std::tr1::bind(&PLYLoader::green_callback, this,  _1);
		} else if (property_name == "diffuse_blue" || property_name == "blue") {
			return std::tr1::bind(&PLYLoader::blue_callback, this,  _1);
		}
	}
	return 0;
}

template<> std::tr1::function <void (ply::uint8)> 
	PLYLoader::scalar_property_definition_callback(const std::string& element_name, 
	const std::string& property_name) {
	if (element_name == "vertex") {
		if (property_name == "diffuse_red" || property_name == "red") {
			return std::tr1::bind(&PLYLoader::red_callback_uint8, this,  _1);
		} else if (property_name == "diffuse_green" || property_name == "green") {
			return std::tr1::bind(&PLYLoader::green_callback_uint8, this,  _1);
		} else if (property_name == "diffuse_blue" || property_name == "blue") {
			return std::tr1::bind(&PLYLoader::blue_callback_uint8, this,  _1);
		}
	}
	return 0;
}

template<> std::tr1::tuple<std::tr1::function<void (ply::uint8)>, 
	std::tr1::function<void (ply::int32)>, std::tr1::function<void ()> > 
	PLYLoader::list_property_definition_callback(const std::string& element_name,
	const std::string& property_name) {
	if ((element_name == "face") && (property_name == "vertex_indices")) {
		return std::tr1::tuple<std::tr1::function<void (ply::uint8)>, 
			std::tr1::function<void (ply::int32)>, std::tr1::function<void ()> >(
			std::tr1::bind(&PLYLoader::face_vertex_indices_begin_uint8, this, _1),
			std::tr1::bind(&PLYLoader::face_vertex_indices_element_int32, this, _1),
			std::tr1::bind(&PLYLoader::face_vertex_indices_end, this)
		);
	} else {
		return std::tr1::tuple<std::tr1::function<void (ply::uint8)>, 
			std::tr1::function<void (ply::int32)>, 
			std::tr1::function<void ()> >(0, 0, 0);
	}
}

template<> std::tr1::tuple<std::tr1::function<void (ply::uint32)>, 
	std::tr1::function<void (ply::int32)>, std::tr1::function<void ()> > 
	PLYLoader::list_property_definition_callback(const std::string& element_name,
	const std::string& property_name) {
	if ((element_name == "face") && (property_name == "vertex_indices")) {
		return std::tr1::tuple<std::tr1::function<void (ply::uint32)>, 
			std::tr1::function<void (ply::int32)>, std::tr1::function<void ()> >(
			std::tr1::bind(&PLYLoader::face_vertex_indices_begin_uint32, this, _1),
			std::tr1::bind(&PLYLoader::face_vertex_indices_element_int32, this, _1),
			std::tr1::bind(&PLYLoader::face_vertex_indices_end, this)
		);
	} else {
		return std::tr1::tuple<std::tr1::function<void (ply::uint32)>, 
			std::tr1::function<void (ply::int32)>, 
			std::tr1::function<void ()> >(0, 0, 0);
	}
}

template<> std::tr1::tuple<std::tr1::function<void (ply::uint8)>, 
	std::tr1::function<void (ply::uint32)>, std::tr1::function<void ()> > 
	PLYLoader::list_property_definition_callback(const std::string& element_name,
	const std::string& property_name) {
	if ((element_name == "face") && (property_name == "vertex_indices")) {
		return std::tr1::tuple<std::tr1::function<void (ply::uint8)>, 
			std::tr1::function<void (ply::uint32)>, std::tr1::function<void ()> >(
			std::tr1::bind(&PLYLoader::face_vertex_indices_begin_uint8, this, _1),
			std::tr1::bind(&PLYLoader::face_vertex_indices_element_uint32, this, _1),
			std::tr1::bind(&PLYLoader::face_vertex_indices_end, this)
		);
	} else {
		return std::tr1::tuple<std::tr1::function<void (ply::uint8)>, 
			std::tr1::function<void (ply::uint32)>, 
			std::tr1::function<void ()> >(0, 0, 0);
	}
}

template<> std::tr1::tuple<std::tr1::function<void (ply::uint32)>, 
	std::tr1::function<void (ply::uint32)>, std::tr1::function<void ()> > 
	PLYLoader::list_property_definition_callback(const std::string& element_name,
	const std::string& property_name) {
	if ((element_name == "face") && (property_name == "vertex_indices")) {
		return std::tr1::tuple<std::tr1::function<void (ply::uint32)>, 
			std::tr1::function<void (ply::uint32)>, std::tr1::function<void ()> >(
			std::tr1::bind(&PLYLoader::face_vertex_indices_begin_uint32, this, _1),
			std::tr1::bind(&PLYLoader::face_vertex_indices_element_uint32, this, _1),
			std::tr1::bind(&PLYLoader::face_vertex_indices_end, this)
		);
	} else {
		return std::tr1::tuple<std::tr1::function<void (ply::uint32)>, 
			std::tr1::function<void (ply::uint32)>, 
			std::tr1::function<void ()> >(0, 0, 0);
	}
}


void PLYLoader::loadPLY(const fs::path &path) {
	ply::ply_parser ply_parser;
	ply_parser.info_callback(std::tr1::bind(&PLYLoader::info_callback,
		this, std::tr1::ref(m_name), _1, _2));
	ply_parser.warning_callback(std::tr1::bind(&PLYLoader::warning_callback, 
		this, std::tr1::ref(m_name), _1, _2));
	ply_parser.error_callback(std::tr1::bind(&PLYLoader::error_callback, 
		this, std::tr1::ref(m_name), _1, _2)); 

	ply_parser.element_definition_callback(std::tr1::bind(&PLYLoader::element_definition_callback,
		this, _1, _2));

	ply::ply_parser::scalar_property_definition_callbacks_type scalar_property_definition_callbacks;
	ply::ply_parser::list_property_definition_callbacks_type list_property_definition_callbacks;

	ply::at<ply::float32>(scalar_property_definition_callbacks) = std::tr1::bind(
		&PLYLoader::scalar_property_definition_callback<ply::float32>, this, _1, _2);

	ply::at<ply::uint8>(scalar_property_definition_callbacks) = std::tr1::bind(
		&PLYLoader::scalar_property_definition_callback<ply::uint8>, this, _1, _2);

	ply::at<ply::uint8, ply::int32>(list_property_definition_callbacks) = std::tr1::bind(
		&PLYLoader::list_property_definition_callback<ply::uint8, ply::int32>, this, _1, _2);

	ply::at<ply::uint32, ply::int32>(list_property_definition_callbacks) = std::tr1::bind(
		&PLYLoader::list_property_definition_callback<ply::uint32, ply::int32>, this, _1, _2);
	
	ply::at<ply::uint8, ply::uint32>(list_property_definition_callbacks) = std::tr1::bind(
		&PLYLoader::list_property_definition_callback<ply::uint8, ply::uint32>, this, _1, _2);

	ply::at<ply::uint32, ply::uint32>(list_property_definition_callbacks) = std::tr1::bind(
		&PLYLoader::list_property_definition_callback<ply::uint32, ply::uint32>, this, _1, _2);

	ply_parser.scalar_property_definition_callbacks(scalar_property_definition_callbacks);
	ply_parser.list_property_definition_callbacks(list_property_definition_callbacks);

	ref<Timer> timer = new Timer();
	ply_parser.parse(path.file_string());

	Log(EInfo, "\"%s\": Loaded " SIZE_T_FMT " triangles, " SIZE_T_FMT 
			" vertices (%s in %i ms).", m_name.c_str(), m_triangleCount, m_vertexCount,
			memString(sizeof(uint32_t) * m_triangleCount * 3 + sizeof(Vertex) * m_vertexCount).c_str(),
			timer->getMilliseconds());
}


MTS_IMPLEMENT_CLASS_S(PLYLoader, false, TriMesh)
MTS_EXPORT_PLUGIN(PLYLoader, "PLY mesh loader");
MTS_NAMESPACE_END
