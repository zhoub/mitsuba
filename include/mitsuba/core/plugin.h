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

#if !defined(__PLUGIN_H)
#define __PLUGIN_H

#include <mitsuba/mitsuba.h>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

MTS_NAMESPACE_BEGIN

/**
 * \brief Abstract plugin class -- represents loadable configurable objects
 * and utilities.
 * 
 * Please see the <tt>\ref ConfigurableObject</tt> and
 * <tt>\ref Utility</tt> classes for details
 */
class MTS_EXPORT_CORE Plugin {
	typedef void *(*CreateInstanceFunc)(const Properties &props);
	typedef void *(*CreateUtilityFunc)();
	typedef char *(*GetDescriptionFunc)();
public:
	/// Load a plugin from the supplied path
	Plugin(const std::string &shortName, const fs::path &path);

	/// Virtual destructor
	virtual ~Plugin();
	
	/// Is this a configurable object plugin or an utility plugin?
	inline bool isUtility() const { return m_isUtility; }

	/// Return an instance of the class implemented by this plugin 
	ConfigurableObject *createInstance(const Properties &props) const;

	/// Return an utility instance (if this is an utility plugin)
	Utility *createUtility() const;

	/// Return a description of this plugin
	std::string getDescription() const;
	
	/// Return the path of this plugin
	inline const fs::path &getPath() const { return m_path; }
	
	/// Return a short name of this plugin
	inline const std::string &getShortName() const { return m_shortName; }
protected:
	/// Resolve the given symbol and return a pointer
	void *getSymbol(const std::string &sym);
	/// Check whether a certain symbol is provided by the plugin
	bool hasSymbol(const std::string &sym) const;
private:
#if defined(WIN32)
	HMODULE m_handle;
#else
	void *m_handle;
#endif
	std::string m_shortName;
	fs::path m_path;
	bool m_isUtility;
	GetDescriptionFunc m_getDescription;
	CreateInstanceFunc m_createInstance;
	CreateUtilityFunc m_createUtility;
};


/**
 * \brief The plugin manager is responsible for resolving and
 * loading external plugins.
 */
class MTS_EXPORT_CORE PluginManager : public Object {
public:
	/// Return the global plugin manager
	inline static PluginManager *getInstance() {
		return m_instance;
	}

	/// Ensure that a plugin is loaded and ready
	void ensurePluginLoaded(const std::string &name);

	/// Return the list of loaded plugins
	std::vector<std::string> getLoadedPlugins() const;

	/**
	 * \brief Instantiate an object using a plugin
	 * \param classType Expected type of the plugin. An
	 *    exception will be thrown if it turns out not
	 *    to derive from this class.
	 * \param props A \ref Properties instance containing
	 *    all information required to find and construct 
	 *    the plugin.
	 */
	ConfigurableObject *createObject(
		const Class *classType,
		const Properties &props
	);

	/// Initializes the global plugin manager instance
	static void staticInitialization();

	/// Free the memory taken by staticInitialization()
	static void staticShutdown();

	MTS_DECLARE_CLASS()
protected:
	PluginManager();
	
	/// Destruct and unload all plugins
	~PluginManager();
private:
	std::map<std::string, Plugin *> m_plugins;
	mutable ref<Mutex> m_mutex;
	static ref<PluginManager> m_instance;
};

MTS_NAMESPACE_END

#endif /* __PLUGIN_H */
