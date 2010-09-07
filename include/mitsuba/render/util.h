#if !defined(__UTILITY_H)
#define __UTILITY_H

#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/** \brief Abstract utility class -- can be used to implement
 * loadable utility plugins that perform various actions. They
 * can be started using the 'mtsutil' launcher.
 */
class MTS_EXPORT_RENDER Utility : public Object {
public:
	/**
	 * Run the utility. The supplied <tt>argc</tt>
	 * and <tt>argv</tt> parameters contain any 
	 * extra arguments passed to mtsutil. The value
	 * returned here will be used as the return value of the
	 * 'mtsutil' process.
	 */
	virtual int run(int argc, char **argv) = 0;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~Utility() { }

	/// Load a scene
	ref<Scene> loadScene(const std::string &fname);
};

#define MTS_DECLARE_UTILITY() \
	MTS_DECLARE_CLASS()

#define MTS_EXPORT_UTILITY(name, descr) \
	MTS_IMPLEMENT_CLASS(name, false, Utility) \
	extern "C" { \
		void MTS_EXPORT *CreateUtility() { \
			return new name(); \
		} \
		const char MTS_EXPORT *GetDescription() { \
			return descr; \
		} \
	}

MTS_NAMESPACE_END

#endif /* __UTILITY_H */
