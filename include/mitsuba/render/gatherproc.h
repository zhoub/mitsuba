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

#if !defined(__GATHERPHOTONPROC_H)
#define __GATHERPHOTONPROC_H

#include <mitsuba/render/particleproc.h>
#include <mitsuba/render/photonmap.h>

MTS_NAMESPACE_BEGIN

/**
 * Process for parallel photon map construction. Given a number and 
 * type (surface/caustic/volume) of photons, it distributes the work
 * over an arbitrary number of machines.
 */
class MTS_EXPORT_RENDER GatherPhotonProcess : public ParticleProcess {
public:
	enum EGatherType {
		/// Surface photons (indirect on diffuse surfaces, last bounce was not through a delta BSDF)
		ESurfacePhotons,
		/// Caustic photons (indirect on diffuse surfaces, last bounce was through a delta BSDF)
		ECausticPhotons,
		/// Surface photons (all of them, even direct)
		EAllSurfacePhotons,
		/// Volumetric photons (multiple scattering only)
		EVolumePhotons
	};

	/**
	 * Create a new process for parallel photon gathering
	 * @param type
	 *     Specifies the type of requested photons (surface/caustic/volume)
	 * @param photonCount
	 *     Specifies the number of requested photons
	 * @param granularity
	 *     Size of the internally used work units (in photons)
	 * @param progressReporterPayload
	 *    Custom pointer payload to be delivered with progress messages
	 */
	GatherPhotonProcess(EGatherType type, size_t photonCount, 
		unsigned int granularity, int maxDepth, int rrDepth,
		const void *progressReporterPayload);

	/**
	 * Once the process has finished, this returns a reference 
	 * to the (still unbalanced) photon map
	 */
	inline PhotonMap *getPhotonMap() { return m_photonMap.get(); }

	/**
	 * Due to asynchronous processing, some excess photons
	 * will generally be produced. This function returns the number
	 * of excess photons that had to be discarded. If this is too
	 * high, the granularity should be decreased.
	 */
	inline size_t getExcess() const { return m_excess; }
	
	/**
	 * Lists the nuber of photons, which had to be shot
	 * in order to fill the photon map.
	 */
	inline size_t getShotPhotons() const { return m_numShot; }

	/* ParallelProcess implementation */
	ref<WorkProcessor> createWorkProcessor() const; 
	void processResult(const WorkResult *wr, bool cancelled);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GatherPhotonProcess() { }
protected:
	EGatherType m_type;
	ref<PhotonMap> m_photonMap;
	int m_maxDepth;
	int m_rrDepth;
	size_t m_excess, m_numShot;
};

MTS_NAMESPACE_END

#endif /* __GATHERPHOTONPROC_H */
