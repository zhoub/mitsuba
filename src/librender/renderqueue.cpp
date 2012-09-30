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

#include <mitsuba/core/timer.h>
#include <mitsuba/render/renderjob.h>

MTS_NAMESPACE_BEGIN

RenderQueue::RenderQueue() {
	m_mutex = new Mutex();
	m_joinMutex = new Mutex();
	m_cond = new ConditionVariable(m_mutex);
	m_timer = new Timer();
}
	
RenderQueue::~RenderQueue() {
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->decRef();
}

void RenderQueue::addJob(RenderJob *job) {
	LockGuard lock(m_mutex);
	m_jobs[job] = JobRecord(m_timer->getMilliseconds());
	job->incRef();
}

void RenderQueue::registerListener(RenderListener *listener) {
	listener->incRef();
	LockGuard lock(m_mutex);
	m_listeners.push_back(listener);
}

void RenderQueue::unregisterListener(RenderListener *listener) {
	{
		LockGuard lock(m_mutex);
		m_listeners.erase(std::remove(m_listeners.begin(), m_listeners.end(),
			listener));
	}
	listener->decRef();
}
	
void RenderQueue::flush() {
	LockGuard lock(m_mutex);
	std::map<RenderJob *, JobRecord>::iterator it = m_jobs.begin();
	for (; it != m_jobs.end(); ++it) {
		(*it).first->flush();
	}
}

void RenderQueue::removeJob(RenderJob *job, bool cancelled) {
	LockGuard lock(m_mutex);
	std::map<RenderJob *, JobRecord>::iterator it = m_jobs.find(job);
	if (it == m_jobs.end()) {
		Log(EError, "RenderQueue::removeRenderJob() - job not found!");
	}
	JobRecord &rec = (*it).second;
	unsigned int ms = m_timer->getMilliseconds() - rec.startTime;
	Log(EInfo, "Render time: %s", timeString(ms/1000.0f, true).c_str());
	m_jobs.erase(job);
	m_cond->broadcast();
	{
		LockGuard lockJoin(m_joinMutex);
		m_joinList.push_back(job);
	}
	signalFinishJob(job, cancelled);
}
	
void RenderQueue::waitLeft(size_t njobs) const {
	UniqueLock lock(m_mutex);
	while (m_jobs.size() > njobs) 
		m_cond->wait();
	lock.unlock();
	join();
}

void RenderQueue::join() const {
	LockGuard lock(m_joinMutex);
	/* Wait for the proper termination of all stopping threads */
	for (size_t i=0; i<m_joinList.size(); ++i) {
		RenderJob *job = m_joinList[i];
		job->join();
		job->decRef();
	}
	m_joinList.clear();
}

void RenderQueue::signalWorkBegin(const RenderJob *job, const RectangularWorkUnit *wu, int worker) {
	LockGuard lock(m_mutex);
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->workBeginEvent(job, wu, worker);
}

void RenderQueue::signalWorkEnd(const RenderJob *job, const ImageBlock *wr) {
	LockGuard lock(m_mutex);
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->workEndEvent(job, wr);
}

void RenderQueue::signalWorkCanceled(const RenderJob *job, const Point2i &offset, const Vector2i &size) {
	LockGuard lock(m_mutex);
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->workCanceledEvent(job, offset, size);
}

void RenderQueue::signalFinishJob(const RenderJob *job, bool cancelled) {
	LockGuard lock(m_mutex);
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->finishJobEvent(job, cancelled);
}

void RenderQueue::signalRefresh(const RenderJob *job) {
	LockGuard lock(m_mutex);
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->refreshEvent(job);
}

MTS_IMPLEMENT_CLASS(RenderQueue, false, Object)
MTS_IMPLEMENT_CLASS(RenderListener, true, Object)
MTS_NAMESPACE_END
