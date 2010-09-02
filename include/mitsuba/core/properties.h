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

#if !defined(__PROPERTIES_H)
#define __PROPERTIES_H

#include <mitsuba/core/transform.h>
#include <mitsuba/core/spectrum.h>
#include <map>

MTS_NAMESPACE_BEGIN

/** \brief Associative map for values of various types. Used to
 * construct sub-classes of <tt>ConfigurableObject</tt>
 */
class MTS_EXPORT_CORE Properties {
public:
	enum PropertyType {
		EBoolean = 0,
		EInteger,
		EFloat,
		EPoint,
		ETransform,
		ESpectrum,
		EString
	};

	/// Construct an empty property container
	Properties() : m_id("unnamed") { }

	/// Construct an empty property container and set the plugin name
	Properties(const std::string &pluginName) : m_pluginName(pluginName), m_id("unnamed") { }

	/// Set the associated plugin name
	inline void setPluginName(const std::string &name) { m_pluginName = name; }
	/// Get the associated plugin name
	inline const std::string &getPluginName() const { return m_pluginName; }
	
	/// Returns the associated ID (or the string "unnamed")
	inline const std::string &getID() const { return m_id; }
	/// Set the associated ID
	inline void setID(const std::string &id) { m_id = id; }

	/// Set a boolean value
	void setBoolean(const std::string &name, bool value, bool warnDuplicates = true);
	/// Get an boolean value
	bool getBoolean(const std::string &name) const;
	/// Get an boolean value (with default);
	bool getBoolean(const std::string &name, bool defVal) const;

	/// Set an integer value
	void setInteger(const std::string &name, int value, bool warnDuplicates = true);
	/// Get an integer value
	int getInteger(const std::string &name) const;
	/// Get an integer value (with default);
	int getInteger(const std::string &name, int defVal) const;

	/// Set an integer value
	void setLong(const std::string &name, int64_t value, bool warnDuplicates = true);
	/// Get an intteger value
	int64_t getLong(const std::string &name) const;
	/// Get an integer value (with default);
	int64_t getLong(const std::string &name, int64_t defVal) const;

	/// Set a single precision floating point value
	void setFloat(const std::string &name, Float value, bool warnDuplicates = true);
	/// Get a single precision floating point value
	Float getFloat(const std::string &name) const;
	/// Get a single precision floating point value (with default)
	Float getFloat(const std::string &name, Float defVal) const;

	/// Set a linear transformation
	void setTransform(const std::string &name, const Transform &value, bool warnDuplicates = true);
	/// Get a linear transformation
	Transform getTransform(const std::string &name) const;
	/// Get a linear transformation (with default)
	Transform getTransform(const std::string &name, const Transform &defVal) const;

	/// Set a spectral power distribution
	void setSpectrum(const std::string &name, const Spectrum &value, bool warnDuplicates = true);
	/// Get a spectral power distribution
	Spectrum getSpectrum(const std::string &name) const;
	/// Get a spectral power distribution (with default)
	Spectrum getSpectrum(const std::string &name, const Spectrum &defVal) const;
	
	/// Set a 3d point
	void setPoint(const std::string &name, const Point &value, bool warnDuplicates = true);
	/// Get a 3d point
	Point getPoint(const std::string &name) const;
	/// Get a 3d point (with default)
	Point getPoint(const std::string &name, const Point &defVal) const;
	/// Get a 3d vector 
	Vector getVector(const std::string &name) const;
	/// Get a 3d vector (with default)
	Vector getVector(const std::string &name, const Vector &defVal) const;

	/// Set a string
	void setString(const std::string &name, const std::string &value, bool warnDuplicates = true);
	/// Get a string
	std::string getString(const std::string &name) const;
	/// Get a string (with default)
	std::string getString(const std::string &name, const std::string &defVal) const;

	/// Store an array containing the names of all stored properties
	inline void putPropertyNames(std::vector<std::string> &results) const {
		for (std::map<std::string, Element>::const_iterator it = m_elements.begin();
			it != m_elements.end(); ++it) 
			results.push_back((*it).first);
	}

	/// Verify if a value with the specified name exists
	bool hasProperty(const std::string &name) const;

	/// Return the property of a type
	PropertyType getType(const std::string &name) const;

	/// Return the list of un-queried attributed
	std::vector<std::string> getUnqueried() const;

	/// Return a string representation
	std::string toString() const;
private:
	struct Element {
		PropertyType type;
		union {
			bool v_boolean;
			int64_t v_long;
			Float v_float;
		};
		// not allowed in union (constructor)
		Point v_point; 
		Vector v_vector;
		Transform v_transform;
		Spectrum v_spectrum;
		std::string v_string;
		mutable bool queried;
	};

	std::map<std::string, Element> m_elements;
	std::string m_pluginName, m_id;
};

MTS_NAMESPACE_END

#endif /* __PROPERTIES_H */
