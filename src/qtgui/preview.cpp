#include "glwidget.h"
#include "preview.h"

PreviewThread::PreviewThread(Device *parentDevice, Renderer *parentRenderer)
	: Thread("prev"), m_parentDevice(parentDevice), m_parentRenderer(parentRenderer), 
		m_context(NULL), m_quit(false) {
	MTS_AUTORELEASE_BEGIN()
	m_session = Session::create();
	m_device = Device::create(m_session);
	m_renderer = Renderer::create(m_session);
	m_mutex = new Mutex();
	m_queueCV = new ConditionVariable(m_mutex);
	m_random = new Random();
	m_bufferCount = 3;
	m_queueEntryIndex = 0;
	m_session->init();
	m_timer = new Timer();
	m_accumBuffer = NULL;
	m_sleep = false;
	m_started = new WaitFlag();
	m_fileResolver = FileResolver::getInstance();

	m_accumProgram = m_renderer->createGPUProgram("Accumulation program");
	m_accumProgram->setSource(GPUProgram::EVertexProgram,
		"void main() {\n"
		"	gl_Position = ftransform();\n"
		"   gl_TexCoord[0]  = gl_MultiTexCoord0;\n"
		"}\n"
	);

	m_accumProgram->setSource(GPUProgram::EFragmentProgram,
		"uniform sampler2D source1, source2;\n"
		"void main() {\n"
		"	gl_FragColor = texture2D(source1, gl_TexCoord[0].xy) + \n"
		"                  texture2D(source2, gl_TexCoord[0].xy);\n"
		"}\n"
	);
				
	m_framebuffer = m_renderer->createGPUTexture("Framebuffer");
	for (int i=0; i<m_bufferCount; ++i) 
		m_recycleQueue.push_back(PreviewQueueEntry(m_queueEntryIndex++));

	MTS_AUTORELEASE_END()
}

PreviewThread::~PreviewThread() {
	MTS_AUTORELEASE_BEGIN()
	m_session->shutdown();
	MTS_AUTORELEASE_END()
}

void PreviewThread::quit() {
	if (!isRunning())
		return;

	std::vector<PreviewQueueEntry> temp;
	temp.reserve(m_bufferCount);
	
	/* Steal all buffers */
	m_mutex->lock();
	while (true) {
		while (m_readyQueue.size() > 0) {
			temp.push_back(m_readyQueue.back());
			m_readyQueue.pop_back();
		}

		while (m_recycleQueue.size() > 0) {
			temp.push_back(m_recycleQueue.back());
			m_recycleQueue.pop_back();
		}

		if ((int) temp.size() == m_bufferCount)
			break;

		m_queueCV->wait();
	}

	/* Put back all buffers and disassociate */
	for (size_t i=0; i<temp.size(); ++i) {
		if (temp[i].buffer)
			temp[i].buffer->disassociate();
		m_recycleQueue.push_back(temp[i]);
	}
	m_quit = true;

	m_queueCV->signal();
	m_mutex->unlock();

	/* Wait for the thread to terminate */
	if (isRunning())
		join();
}

void PreviewThread::setSceneContext(SceneContext *context, bool swapContext, bool motion) {
	std::vector<PreviewQueueEntry> temp;
	temp.reserve(m_bufferCount);

	m_sleep = true;
	m_mutex->lock();

	/* Steal all buffers from the rendering
	   thread to make sure we get its attention :) */
	while (true) {
		while (m_readyQueue.size() > 0) {
			temp.push_back(m_readyQueue.front());
				m_readyQueue.pop_front();
		}

		while (m_recycleQueue.size() > 0) {
			temp.push_back(m_recycleQueue.back());
			m_recycleQueue.pop_back();
		}

		if ((int) temp.size() == m_bufferCount)
			break;

		m_queueCV->wait();
	}

	if (swapContext && m_context) {
		m_context->vpls = m_vpls;
		m_context->previewBuffer = temp[0];
		m_context->previewBuffer.buffer->disassociate();
		m_recycleQueue.push_back(PreviewQueueEntry(m_queueEntryIndex++));
		
		/* Put back all buffers */
		for (size_t i=1; i<temp.size(); ++i)
			m_recycleQueue.push_back(temp[i]);
	} else {
		for (size_t i=0; i<temp.size(); ++i)
			m_recycleQueue.push_back(temp[i]);
	}

	if (swapContext && context && context->previewBuffer.vplSampleOffset > 0) {
		/* Resume from a stored state */
		m_vplSampleOffset = context->previewBuffer.vplSampleOffset;
		m_vpls = context->vpls;
		m_accumBuffer = context->previewBuffer.buffer;

		/* Take ownership of the buffer */
		m_recycleQueue.push_back(context->previewBuffer);
		context->previewBuffer.buffer = NULL;
		context->previewBuffer.sync = NULL;
		context->previewBuffer.vplSampleOffset = 0;

		if (m_recycleQueue.size() > (size_t) m_bufferCount) {
			PreviewQueueEntry entry = m_recycleQueue.front();
			m_recycleQueue.pop_front();
			if (entry.buffer) {
				entry.buffer->disassociate();
				entry.buffer->decRef();
			}
			if (entry.sync) 
				entry.sync->decRef();
		}
	} else {
		/* Reset the VPL rendering progress */
		m_vplSampleOffset = 0;
		m_vpls.clear();
		m_accumBuffer = NULL;
	}

	if (m_context != context)
		m_minVPLs = 0;

	m_vplsPerSecond = 0;
	m_raysPerSecond = 0;
	m_vplCount = 0;
	m_timer->reset();	

	m_context = context;

	if (m_context) {
		Camera *camera = m_context->scene->getCamera();
		m_camPos = camera->getPosition();
		m_camViewTransform = camera->getViewTransform();
	}

	if (motion && !m_motion) {
		emit statusMessage("");
		m_minVPLs = 0;
	}

	m_motion = motion;
	m_queueCV->signal();
	m_mutex->unlock();
	m_sleep = false;
}

void PreviewThread::resume() {
	m_queueCV->signal();
}

PreviewQueueEntry PreviewThread::acquireBuffer(int ms) {
	PreviewQueueEntry entry;

	m_mutex->lock();
	while (m_readyQueue.size() == 0) {
		if (m_quit)
			return entry;
		if (!m_queueCV->wait(ms)) {
			m_mutex->unlock();
			return entry;
		}
	}
	entry = m_readyQueue.front();
	m_readyQueue.pop_front();
	m_mutex->unlock();

	if (m_context->previewMethod == ERayTrace || 
		m_context->previewMethod == ERayTraceCoherent) 
		entry.buffer->refresh();
	else if (m_useSync) 
		entry.sync->enqueueWait(); 

	return entry;
}

void PreviewThread::releaseBuffer(PreviewQueueEntry &entry) {
	m_mutex->lock();

	if (m_motion)
		m_readyQueue.push_front(entry);
	else
		m_recycleQueue.push_back(entry);

	if (m_useSync)
		entry.sync->cleanup();

	m_queueCV->signal();
	m_mutex->unlock();
}

void PreviewThread::run() {
	MTS_AUTORELEASE_BEGIN()

	bool initializedGraphics = false;
	FileResolver::setInstance(m_fileResolver);

	try {
		m_device->init(m_parentDevice);
		m_device->setVisible(false);

		/* We have alrady seen this once */
		m_renderer->setLogLevel(ETrace);
		m_renderer->setWarnLogLevel(ETrace);
		m_renderer->init(m_device, m_parentRenderer);
		m_renderer->setLogLevel(EDebug);
		m_renderer->setWarnLogLevel(EWarn);
		m_started->set(true);

		m_accumProgram->init();
		m_accumProgramParam_source1 = m_accumProgram->getParameterID("source1");
		m_accumProgramParam_source2 = m_accumProgram->getParameterID("source2");
		m_useSync = m_renderer->getCapabilities()->isSupported(RendererCapabilities::ESyncObjects);

		initializedGraphics = true;

		while (true) {
			PreviewQueueEntry target;

			m_mutex->lock();
			while (!(m_quit || (m_context != NULL && m_context->mode == EPreview 
					&& ((m_readyQueue.size() != 0 && !m_motion) || m_recycleQueue.size() != 0))))
				m_queueCV->wait();

			MTS_AUTORELEASE_END()
			MTS_AUTORELEASE_BEGIN()

			if (m_quit) {
				m_mutex->unlock();
				break;
			} else if (m_recycleQueue.size() != 0) {
				target = m_recycleQueue.front();
				m_recycleQueue.pop_front();
			} else if (m_readyQueue.size() != 0 && !m_motion) {
				target = m_readyQueue.front();
				m_readyQueue.pop_front();
			} else {
				Log(EError, "Internal error!");
			}

			if (m_motion && m_vplCount >= m_minVPLs && m_minVPLs != 0) {
				/* The user is currently moving around, and a good enough
				   preview has already been rendered. Don't improve it to
				   avoid flicker */
				m_recycleQueue.push_back(target);
				m_queueCV->wait();
				m_mutex->unlock();
				continue;
			}

			m_mutex->unlock();

			if (m_vplSampleOffset == 0) 
				m_accumBuffer = NULL;

			const Film *film = m_context->scene->getFilm();
			Point3i size(film->getSize().x, film->getSize().y, 1);

			if (target.buffer == NULL || target.buffer->getSize() != size) {
				target.buffer = m_renderer->createGPUTexture(formatString("Communication buffer %i", target.id));
				target.buffer->setFormat(GPUTexture::EFloat32RGB);
				target.buffer->setSize(size);
				target.buffer->setFilterType(GPUTexture::ENearest);
				target.buffer->setFrameBufferType(GPUTexture::EColorBuffer);
				target.buffer->setMipMapped(false);
				target.buffer->init();
				target.buffer->incRef();
				target.sync = m_renderer->createGPUSync();
				target.sync->incRef();
			}

			if (m_context->previewMethod == ERayTrace || m_context->previewMethod == ERayTraceCoherent) {
				if (m_previewProc == NULL || m_previewProc->getScene() != m_context->scene) 
					m_previewProc = new PreviewProcess(m_context->scene, m_context->sceneResID, 32);

				if (m_timer->getMilliseconds() > 1000) {
					Float time = 1000 / (Float) m_timer->getMilliseconds();
					Float vplCount = m_vplsPerSecond * time, rayCount = m_raysPerSecond * time / 1000000.0f;
					if (!m_motion)
						emit statusMessage(QString(formatString("%.1f VPLs, %.1f MRays/sec", vplCount, rayCount).c_str()));
					m_vplsPerSecond = 0;
					m_raysPerSecond = 0;
					m_timer->reset();
				}

				if (m_vpls.empty())
					m_vplSampleOffset = generateVPLs(m_context->scene, m_vplSampleOffset,
						1, m_context->pathLength, m_vpls);

				VPL vpl = m_vpls.front();
				m_vpls.pop_front();

				rtrtRenderVPL(target, vpl);
			} else {
				if (m_shaderManager == NULL || m_shaderManager->getScene() != m_context->scene) {
					if (m_shaderManager) {
						m_shaderManager->cleanup();
						m_framebuffer->cleanup();
					}
					m_shaderManager = new VPLShaderManager(m_context->scene, m_renderer);
					m_shaderManager->init();
					m_framebuffer->setFormat(GPUTexture::EFloat32RGB);
					m_framebuffer->setSize(size);
					m_framebuffer->setFilterType(GPUTexture::ENearest);
					m_framebuffer->setFrameBufferType(GPUTexture::EColorBuffer);
					m_framebuffer->setMipMapped(false);
					m_framebuffer->init();
				}

				m_shaderManager->setShadowMapResolution(m_context->shadowMapResolution);
				m_shaderManager->setClamping(m_context->clamping);
				m_shaderManager->setSinglePass(m_context->previewMethod == EOpenGLSinglePass);

				if (m_timer->getMilliseconds() > 1000) {
					Float count = m_vplsPerSecond / (Float) m_timer->getMilliseconds() * 1000;
					if (!m_motion)
						emit statusMessage(QString(formatString("%.1f VPLs/sec", count).c_str()));
					m_vplsPerSecond = 0;
					m_timer->reset();
				}

				if (m_vpls.empty())
					m_vplSampleOffset = generateVPLs(m_context->scene, m_vplSampleOffset,
						1, m_context->pathLength, m_vpls);

				VPL vpl = m_vpls.front();
				m_vpls.pop_front();

				oglRenderVPL(target, vpl);
				
				if (m_useSync)
					target.sync->init(); 
			}

			m_mutex->lock();
			m_vplsPerSecond++;
			m_vplCount++;

			if (m_minVPLs == 0) {
				if (m_timer->getMilliseconds() > 50) 
					m_minVPLs = m_vplCount;
			}

			if (m_vplCount >= m_minVPLs && m_minVPLs > 0)
				m_readyQueue.push_back(target);
			else
				m_recycleQueue.push_back(target);
			m_queueCV->signal();
			m_mutex->unlock();

			if (m_sleep)
				sleep(10);
		}
	} catch (std::exception &e) {
		m_started->set(true);
		Log(EWarn, "Caught an exception: %s", e.what());
		emit caughtException(e.what());
	}

	if (initializedGraphics) {
		if (m_shaderManager)
			m_shaderManager->cleanup();

		m_accumProgram->cleanup();

		m_mutex->lock();
		while (!m_readyQueue.empty()) {
			PreviewQueueEntry &entry = m_readyQueue.back();
			if (entry.buffer) {
				entry.buffer->disassociate();
				entry.buffer->decRef();
			}
			if (entry.sync) 
				entry.sync->decRef();
			m_readyQueue.pop_back();
		}

		while (!m_recycleQueue.empty()) {
			PreviewQueueEntry &entry = m_recycleQueue.back();
			if (entry.buffer) {
				entry.buffer->disassociate();
				entry.buffer->decRef();
			}
			if (entry.sync) 
				entry.sync->decRef();
			m_recycleQueue.pop_back();
		}

		m_renderer->shutdown();
		m_device->shutdown();
		m_mutex->unlock();
	}

	MTS_AUTORELEASE_END()
}

void PreviewThread::oglRenderVPL(PreviewQueueEntry &target, const VPL &vpl) {
	const std::vector<TriMesh *> meshes = m_context->scene->getMeshes();

	m_shaderManager->setVPL(vpl);

	Point2 jitter(.5f, .5f);
	if (!m_motion)
		jitter -= Point2(m_random->nextFloat(), m_random->nextFloat());

	m_mutex->lock();
	const ProjectiveCamera *camera = static_cast<const ProjectiveCamera *>
		(m_context->scene->getCamera());
	Transform projectionTransform = camera->getGLProjectionTransform(jitter);
	m_renderer->setCamera(projectionTransform.getMatrix(), m_camViewTransform.getMatrix());
	Transform clipToWorld = m_camViewTransform.inverse() 
		* Transform::scale(Vector(1, 1, -1)) * projectionTransform.inverse();

	target.vplSampleOffset = m_vplSampleOffset;
	Point camPos = m_camPos;
	m_mutex->unlock();

	m_framebuffer->activateTarget();
	m_framebuffer->clear();
	m_renderer->beginDrawingMeshes();
	for (unsigned int j=0; j<meshes.size(); j++) {
		m_shaderManager->configure(vpl, meshes[j]->getBSDF(), 
			meshes[j]->getLuminaire(), camPos);
		m_renderer->drawTriMesh(meshes[j]);
		m_shaderManager->unbind();
	}
	m_renderer->endDrawingMeshes();
	m_shaderManager->drawBackground(clipToWorld, camPos);
	m_framebuffer->releaseTarget();

	target.buffer->activateTarget();
	m_renderer->setDepthMask(false);
	m_renderer->setDepthTest(false);
	m_framebuffer->bind(0);
	if (m_accumBuffer == NULL) { 
		target.buffer->clear();
		m_renderer->blitTexture(m_framebuffer, true);
	} else {
		m_accumBuffer->bind(1);
		m_accumProgram->bind();
		m_accumProgram->setParameter(m_accumProgramParam_source1, m_accumBuffer);
		m_accumProgram->setParameter(m_accumProgramParam_source2, m_framebuffer);
		m_renderer->blitQuad(true);
		m_accumProgram->unbind();
		m_accumBuffer->unbind();
	}
	m_framebuffer->unbind();
	m_renderer->setDepthMask(true);
	m_renderer->setDepthTest(true);
	target.buffer->releaseTarget();
	m_accumBuffer = target.buffer;

	static int i=0;
	if ((++i % 4) == 0 || m_motion) {
		/* Don't let the queue get too large -- this makes
		   the whole system unresponsive */
		m_renderer->finish();
	} else {
		if (m_useSync) {
			m_renderer->flush();
		} else {
			/* No sync objects available - we have to wait 
			   for everything to finish */
			m_renderer->finish();
		}
	}
}

void PreviewThread::rtrtRenderVPL(PreviewQueueEntry &target, const VPL &vpl) {
	Float nearClip =  std::numeric_limits<Float>::infinity(),
		  farClip  = -std::numeric_limits<Float>::infinity();

	const int sampleCount = 200;
	const Float invSampleCount = 1.0f/sampleCount;
	Ray ray;
	ray.o = vpl.its.p;
	Intersection its;

	for (int i=1; i<=sampleCount; ++i) {
		Vector dir;
		Point2 seed(i*invSampleCount, radicalInverse(2, i)); // Hammersley seq.
		if (vpl.type == ELuminaireVPL && vpl.luminaire->getType() == Luminaire::EDeltaPosition)
			dir = squareToSphere(seed);
		else
			dir = vpl.its.shFrame.toWorld(squareToHemispherePSA(seed));
		ray.setDirection(dir);

		if (m_context->scene->rayIntersect(ray, its)) {
			nearClip = std::min(nearClip, its.t);
			farClip = std::max(farClip, its.t);
		}
	}

	Float minDist = nearClip + (farClip - nearClip) * m_context->clamping;

	if (nearClip >= farClip) {
		/* Unable to find any surface - just default values based on the scene size */
		nearClip = 1e-3 * m_context->scene->getBSphere().radius;
		farClip = 1e3 * m_context->scene->getBSphere().radius;
		minDist = 0;
	}

	Point2 jitter(.5f, .5f);
	if (!m_motion)
		jitter = Point2(m_random->nextFloat(), m_random->nextFloat());

	if (target.buffer->getBitmap() == NULL) 
		target.buffer->setBitmap(0, new Bitmap(target.buffer->getSize().x, 
			target.buffer->getSize().y, 96));

	m_mutex->lock();
	m_previewProc->configure(vpl, minDist, jitter, 
		m_accumBuffer ? m_accumBuffer->getBitmap() : NULL, 
		target.buffer->getBitmap(), m_context->previewMethod == ERayTraceCoherent);
	m_mutex->unlock();

	ref<Scheduler> sched = Scheduler::getInstance();
	sched->schedule(m_previewProc);
	sched->wait(m_previewProc);
	target.vplSampleOffset = m_vplSampleOffset;
	m_raysPerSecond += m_previewProc->getRayCount();
	m_accumBuffer = target.buffer;
}
