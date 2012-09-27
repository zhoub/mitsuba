/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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

#include <mitsuba/bidir/util.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/plugin.h>
#include "pssmlt_proc.h"
#include "pssmlt_sampler.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{pssmlt}{Primary Sample Space Metropolis Light Transport}
 * \order{9}
 * \parameters{
 *	   \parameter{bidirectional}{\Boolean}{
 *	   PSSMLT works in conjunction with another rendering
 *	   technique that is endowed with Markov Chain-based sample generation.
 *	   Two choices are available:
 *	    \begin{itemize}
 *	    \item \code{true}: Operate on top of a fully-fleged bidirectional
 *	      path tracer with multiple importance sampling.
 *	    \item \code{false}: Rely on a unidirectional
 *	    volumetric path tracer (i.e. \pluginref{volpath})
 *	    \vspace{-4mm}
 *	    \end{itemize}
 *	   }
 *     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *         in the generated output image (where \code{-1} corresponds to $\infty$).
 *	       A value of \code{1} will only render directly visible light sources.
 *	       \code{2} will lead to single-bounce (direct-only) illumination, 
 *	       and so on. \default{\code{-1}}
 *	   }
 *	   \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after 
 *	      which the implementation will start to use the ``russian roulette'' 
 *	      path termination criterion. \default{\code{5}}
 *	   }
 *	   \parameter{luminanceSamples}{\Integer}{
 *	      MLT-type algorithms create output images that are only
 *	      \emph{relative}. The algorithm can e.g. determine that a certain pixel
 *	      is approximately twice as bright as another one, but the absolute
 *	      scale is unknown. To recover it, this plugin computes
 *	      the average luminance arriving at the sensor by generating a
 *	      number of samples. \default{\code{100000} samples}
 *     }
 *	   \parameter{pLarge}{\Float}{
 *	     Rate at which the implementation tries to replace the current path
 *	     with a completely new one. Usually, there is little need to change
 *	     this. \default{0.3}
 *	   }
 * }
 * Primary Sample Space Metropolis Light Transport (PSSMLT) is a rendering
 * technique based on Markov Chain Monte Carlo (MCMC) integration. It
 * is based on work by Kelemen et al. \cite{Kelemen2002Simple}.
 *
 * In contrast to simple methods like \pluginref{path} tracing that render
 * images by performing a na\"ive and memoryless random search for light paths, 
 * PSSMLT actively searches for \emph{relevant} light paths (as is the case 
 * for other MCMC methods). Once such a path is found, the algorithm tries to 
 * explore neighboring paths to amortize the cost of the search. This can
 * significantly improve the convergence rate of difficult input.
 * Scenes that were already relatively easy to render usually don't benefit
 * much from PSSMLT, since the MCMC data management causes additional 
 * computational overheads.
 *
 * An interesting aspect of PSSMLT is that it performs this exploration
 * of light paths by perturbing the ``random numbers'' that were initially
 * used to construct the path. After regenerating the path with the 
 * perturbed numbers, it ends up in a slightly different configuration, and
 * this process repeats over and over again.
 * The path regeneration step is fairly arbitrary and this is what makes
 * the method powerful: in particular, it is possible to use PSSMLT as a
 * layer on top of an existing method, which creates a new ``metropolized''
 * rendering algorithm that is enhanced with a certain degree of adaptiveness as 
 * described above.
 *
 * The PSSMLT implementation in Mitsuba can operate on top of either a simple 
 * volumetric path tracer or a fully-fledged bidirectional path tracer with 
 * multiple importance sampling, and this choice is controlled by the 
 * \code{bidirectional} flag.
 *
 * \remarks{
 *    \item This integrator does not work with dipole-style subsurface scattering models.
 * }
 */

class PSSMLT : public Integrator {
public:
	PSSMLT(const Properties &props) : Integrator(props) {
		/* Longest visualized path length (<tt>-1</tt>=infinite). 
		   A value of <tt>1</tt> will visualize only directly visible light 
		   sources. <tt>2</tt> will lead to single-bounce (direct-only) 
		   illumination, and so on. */
		m_config.maxDepth = props.getInteger("maxDepth", -1);

		/* Depth to begin using russian roulette (set to -1 to disable) */
		m_config.rrDepth = props.getInteger("rrDepth", 5);

		/* If set to <tt>true</tt>, the MLT algorithm runs on top of a
		   bidirectional path tracer with multiple importance sampling.
		   Otherwise, the implementation reverts to a basic path tracer.
		   Generally, the bidirectinal path tracer should be noticably
		   better, so it's best to this setting at its default. */
		m_config.technique = props.getBoolean("bidirectional", true) ?		
			PathSampler::EBidirectional : PathSampler::EUnidirectional;

		/* This setting can be very useful to reduce noise in dark regions
		   of the image: it activates two-stage MLT, where a nested MLT renderer
		   first creates a tiny version of the output image. In a second pass, 
		   the full version is then rendered, while making use of information 
		   about the image-space luminance distribution found in the first
		   pass. Two-stage MLT is very useful in making the noise characteristics
		   more uniform over time image -- specifically, since MLT tends to get
		   stuck in very bright regions at the cost of the remainder of the image.*/
		m_config.twoStage = props.getBoolean("twoStage", false);

		/* When running two-stage MLT, this parameter determines the size
		   of the downsampled image created in the first pass (i.e. setting this
		   to 16 means that the horizontal/vertical resolution will be 16 times
		   lower). Usually, it's fine to leave this parameter unchanged. When 
		   the two-stage process introduces noisy halos around very bright image 
		   regions, it can be set to a lower value */
		m_config.firstStageSizeReduction = props.getInteger(
			"firstStageSizeReduction", 16);

		/* Used internally to let the nested rendering process of a
		   two-stage MLT approach know that it is running the first stage */
		m_config.firstStage= props.getBoolean("firstStage", false);

		/* Number of samples used to estimate the total luminance 
		   received by the sensor's sensor */
		m_config.luminanceSamples = props.getInteger("luminanceSamples", 100000);

		/* Probability of creating large mutations in the [Kelemen et. al] 
		   MLT variant. The default is 0.3. */
		m_config.pLarge = props.getFloat("pLarge", 0.3f);

		/* Should direct illumination be handled separately? (i.e. not
		   using MLT) This is usually the right way to go, since direct
		   illumination is easily handled using more optimized rendering
		   techniques that can make use of low-discrepancy point sets.
		   This in turn lets MLT focus on the more difficult parts of the
		   light transport. On the other hand, some scenes use very
		   hard to find paths even for direct illumination, in which case
		   it may make more sense to set this property to 'false' */
		m_config.separateDirect = props.getBoolean("separateDirect",
				true);

		/* When 'separateDirect' is set to 'true', this parameter can
		   be used to specify the samples per pixel used to render the
		   direct component. Should be a power of two (otherwise, it will
		   be rounded to the next one). When set to zero or less, the
		   direct illumination component will be hidden, which is useful
		   for analyzing the component rendered by MLT. */
		m_config.directSamples = props.getInteger("directSamples", 16);

		/* Should the multiple importance sampling-based weight computation by 
		   Kelemen et al. be used? Otherwise, the implementation falls back
		   to the 'use of expectations' technique from Veach-style MLT. */
		m_config.kelemenStyleWeights = props.getBoolean("kelemenStyleWeights", true);

		/* Should an optimized direct illumination sampling strategy be used
		   for s=1 paths? (as opposed to plain emission sampling). Usually
		   a good idea. Note that this setting only applies when the 
		   bidirectional path tracer is used internally. The optimization
		   affects all paths, not just the ones contributing direct illumination,
		   hence it is completely independent of the <tt>separateDirect</tt>
		   parameter. */
		m_config.directSampling = props.getBoolean(
				"directSampling", true);

		/* Recommended mutation sizes in the primary sample space */
		m_config.mutationSizeLow  = props.getFloat("mutationSizeLow",  1.0f/1024.0f);
		m_config.mutationSizeHigh = props.getFloat("mutationSizeHigh", 1.0f/64.0f);
		Assert(m_config.mutationSizeLow > 0 && m_config.mutationSizeHigh > 0 &&
		       m_config.mutationSizeLow < 1 && m_config.mutationSizeHigh < 1 &&
			   m_config.mutationSizeLow < m_config.mutationSizeHigh);

		/* Specifies the number of parallel work units required for
		   multithreaded and network rendering. When set to <tt>-1</tt>, the 
		   amount will default to four times the number of cores. Note that
		   every additional work unit entails a significant amount of 
		   communication overhead (a full-sized floating put image must be 
		   transmitted), hence it is important to set this value as low as 
		   possible, while ensuring that there are enough units to keep all 
		   workers busy. */
		m_config.workUnits = props.getInteger("workUnits", -1);
		m_config.timeout = props.getInteger("timeout", 0);
	}

	/// Unserialize from a binary data stream
	PSSMLT(Stream *stream, InstanceManager *manager)
	 : Integrator(stream, manager) {
		m_config = PSSMLTConfiguration(stream);
		configure();
	}

	virtual ~PSSMLT() { }

	void serialize(Stream *stream, InstanceManager *manager) const {
		Integrator::serialize(stream, manager);
		m_config.serialize(stream);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, 
			const RenderJob *job, int sceneResID, int sensorResID,
			int samplerResID) {
		Integrator::preprocess(scene, queue, job, sceneResID,
				sensorResID, samplerResID);

		ref<const Sensor> sensor = scene->getSensor();

		if (scene->getSubsurfaceIntegrators().size() > 0)
			Log(EError, "Subsurface integrators are not supported by MLT!");

		if (sensor->getSampler()->getClass()->getName() != "IndependentSampler")
			Log(EError, "Metropolis light transport requires the independent sampler");

		return true;
	}

	void cancel() {
		ref<RenderJob> nested = m_nestedJob;
		if (nested)
			nested->cancel();
		Scheduler::getInstance()->cancel(m_process);
	}

	bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
			int sceneResID, int sensorResID, int samplerResID) {
		ref<Scheduler> scheduler = Scheduler::getInstance();
		ref<Sensor> sensor = scene->getSensor();
		ref<Sampler> sampler = sensor->getSampler();
		const Film *film = sensor->getFilm();
		size_t nCores = scheduler->getCoreCount();
		size_t sampleCount = sampler->getSampleCount();
		m_config.importanceMap = NULL;

		if (m_config.twoStage && !m_config.firstStage) {
			Log(EInfo, "Executing first MLT stage");
			ref<Timer> timer = new Timer();
			Assert(m_config.firstStageSizeReduction > 0);
			m_config.importanceMap = BidirectionalUtils::mltLuminancePass(
					scene, sceneResID, queue, m_config.firstStageSizeReduction,
					m_nestedJob);
			if (!m_config.importanceMap) {
				Log(EWarn, "First-stage MLT process failed!");
				return false;
			}
			Log(EInfo, "First MLT stage took %i ms", timer->getMilliseconds());

			std::string debugFile = "mlt_stage1.exr";
			Log(EInfo, "Writing upsampled luminances to \"%s\"", debugFile.c_str());
			ref<FileStream> fs = new FileStream(debugFile, FileStream::ETruncReadWrite);
			m_config.importanceMap->write(Bitmap::EOpenEXR, fs);
		}

		bool nested = m_config.twoStage && m_config.firstStage;

		Vector2i cropSize = film->getCropSize();;
		Log(EInfo, "Starting %srender job (%ix%i, " SIZE_T_FMT
			" %s, " SSE_STR ", approx. " SIZE_T_FMT " mutations/pixel) ..", 
			nested ? "nested " : "", cropSize.x, cropSize.y,
			nCores, nCores == 1 ? "core" : "cores", sampleCount);

		if (m_config.workUnits <= 0) 
			m_config.workUnits = (size_t) std::ceil((cropSize.x * cropSize.y
				* sampleCount) / 200000.0f);

		m_config.nMutations = (cropSize.x * cropSize.y *
			sampleCount) / m_config.workUnits;

		ref<Bitmap> directImage;
		if (m_config.separateDirect && m_config.directSamples > 0 && !nested) {
			directImage = BidirectionalUtils::renderDirectComponent(scene, 
				sceneResID, sensorResID, queue, job, m_config.directSamples);
			if (directImage == NULL)
				return false;
		}

		std::vector<PathSeed> pathSeeds;
		ref<ReplayableSampler> rplSampler = new ReplayableSampler();
		ref<PathSampler> pathSampler = new PathSampler(m_config.technique, scene, 
			rplSampler, rplSampler, rplSampler, m_config.maxDepth, m_config.rrDepth,
			m_config.separateDirect, m_config.directSampling);

		ref<PSSMLTProcess> process = new PSSMLTProcess(job, queue, 
				m_config, directImage, pathSeeds);
	
		m_config.luminance = pathSampler->generateSeeds(m_config.luminanceSamples, 
			m_config.workUnits, false, pathSeeds);

		if (!nested)
			m_config.dump();
		
		/* Create a sampler instance for each worker */
		ref<PSSMLTSampler> mltSampler = new PSSMLTSampler(m_config);
		std::vector<SerializableObject *> mltSamplers(scheduler->getCoreCount());
		for (size_t i=0; i<mltSamplers.size(); ++i) {
			ref<Sampler> clonedSampler = mltSampler->clone();
			clonedSampler->incRef();
			mltSamplers[i] = clonedSampler.get();
		}
		int mltSamplerResID = scheduler->registerMultiResource(mltSamplers);
		for (size_t i=0; i<scheduler->getCoreCount(); ++i)
			mltSamplers[i]->decRef();
		int rplSamplerResID = scheduler->registerResource(rplSampler);

		process->bindResource("scene", sceneResID);
		process->bindResource("sensor", sensorResID);
		process->bindResource("sampler", mltSamplerResID);
		process->bindResource("rplSampler", rplSamplerResID);

		m_process = process;
		scheduler->schedule(process);
		scheduler->wait(process);
		m_process = NULL;
		scheduler->unregisterResource(rplSamplerResID);
		process->develop();
		
		return process->getReturnStatus() == ParallelProcess::ESuccess;
	}

	MTS_DECLARE_CLASS()
private:
	ref<ParallelProcess> m_process;
	ref<RenderJob> m_nestedJob;
	PSSMLTConfiguration m_config;
};

MTS_IMPLEMENT_CLASS_S(PSSMLT, false, Integrator)
MTS_EXPORT_PLUGIN(PSSMLT, "Primary Sample Space MLT");
MTS_NAMESPACE_END
