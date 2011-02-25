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

#if !defined(__CORE_FWD_H)
#define __CORE_FWD_H

MTS_NAMESPACE_BEGIN

class Appender;
class Bitmap;
class BlackBodySpectrum;
template <typename _PointType> struct BoundingBox;
struct BoundingBox3;
struct BoundingSphere;
class ConfigurableObject;
struct CacheLineCounter;
class Class;
class ConditionVariable;
class ConsoleStream;
class DefaultFormatter;
struct DiscretePDF;
class FileResolver;
class FileStream;
class Formatter;
struct Frame;
class InstanceManager;
class InterpolatedSpectrum;
class LocalWorker;
class Logger;
class MemoryStream;
class Mutex;
class NetworkedObject;
struct Normal;
class Object;
class ParallelProcess;
class Plugin;
class PluginManager;
class ProgressReporter;
class Properties;
class Random;
struct Ray;
struct RayDifferential;
class RemoteProcess;
class RemoteWorker;
class RemoteWorkerReader;
class Scheduler;
class Serializable;
class SerializableObject;
class SHRotation;
class SHSampler;
struct SHVector;
struct SHVector4D;
class SmoothSpectrum;
class SocketStream;
class SparseWavelet2D;
class SparseWaveletOctree;
struct Spectrum;
class SSHStream;
class Statistics;
class StatsCounter;
class StreamAppender;
class StreamBackend;
class Thread;
class Timer;
struct Transform;
struct Triangle;
class UnbufferedAppender;
class WaitFlag;
class Wavelet2D;
class Wavelet3D;
class Worker;
class WorkProcessor;
class WorkResult;
class WorkUnit;
struct Quaternion;
class ZStream;

MTS_NAMESPACE_END

#endif /* __CORE_FWD_H */
