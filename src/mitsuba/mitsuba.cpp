/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#include <mitsuba/core/platform.h>
#include <xercesc/parsers/SAXParser.hpp>
#include <mitsuba/core/sched_remote.h>
#include <mitsuba/core/sstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/appender.h>
#include <mitsuba/core/sshstream.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/render/scenehandler.h>
#include <fstream>
#include <stdexcept>

#if defined(WIN32)
#include <mitsuba/core/getopt.h>
#else
#include <signal.h>
#endif

using namespace mitsuba;

void help() {
	cout <<  "Mitsuba version " MTS_VERSION ", Copyright (c) " MTS_YEAR " Wenzel Jakob" << endl;
	cout <<  "Usage: mitsuba [options] <One or more scene XML files>" << endl;
	cout <<  "Options/Arguments:" << endl;
	cout <<  "   -h          Display this help text" << endl << endl;
	cout <<  "   -D key=val  Define a constant, which can referenced as \"$key\" in the scene" << endl << endl;
	cout <<  "   -o fname    Write the output image to the file denoted by \"fname\"" << endl << endl;
	cout <<  "   -a p1;p2;.. Add one or more entries to the resource search path" << endl << endl;
	cout <<  "   -p count    Override the detected number of processors. Useful for reducing" << endl;
	cout <<  "               the load or creating scheduling-only nodes in conjunction with"  << endl;
	cout <<  "               the -c and -s parameters, e.g. -p 0 -c host1;host2;host3,..." << endl << endl;
	cout <<  "   -q          Quiet mode - do not print any log messages to stdout" << endl << endl;
	cout <<  "   -c hosts    Network rendering: connect to mtssrv instances over a network." << endl;
	cout <<  "               Requires a semicolon-separated list of host names of the form" << endl;
	cout <<  "                       host.domain[:port] for a direct connection" << endl;
	cout <<  "                 or" << endl;
	cout <<  "                       user@host.domain[:path] for a SSH connection (where" << endl;
	cout <<  "                       \"path\" denotes the place where Mitsuba is checked" << endl;
	cout <<  "                       out -- by default, \"~/mitsuba\" is used)" << endl << endl;
	cout <<  "   -s file     Connect to additional Mitsuba servers specified in a file" << endl;
	cout <<  "               with one name per line (same format as in -c)" << endl<< endl;
	cout <<  "   -j count    Simultaneously schedule several scenes. Can sometimes accelerate" << endl;
	cout <<  "               rendering when large amounts of processing power are available" << endl;
	cout <<  "               (e.g. when running Mitsuba on a cluster. Default: 1)" << endl << endl;
	cout <<  "   -n name     Assign a node name to this instance (Default: host name)" << endl << endl;
	cout <<  "   -t          Test case mode (see Mitsuba docs for more information)" << endl << endl;
	cout <<  "   -x          Skip rendering of files where output already exists" << endl << endl;
	cout <<  "   -r sec      Write (partial) output images every 'sec' seconds" << endl << endl;
	cout <<  "   -b res      Specify the block resolution used to split images into parallel" << endl;
	cout <<  "               workloads (default: 32). Only applies to some integrators." << endl << endl;
	cout <<  "   -v          Be more verbose" << endl << endl;
	cout <<  "   -w          Treat warnings as errors" << endl << endl;
	cout <<  "   -z          Disable progress bars" << endl << endl;
	cout <<  " The README file included with the distribution contains further information." << endl;
}

ref<RenderQueue> renderQueue = NULL;

#if !defined(WIN32)
/* Handle the hang-up signal and write a partially rendered image to disk */
void signalHandler(int signal) {
	if (signal == SIGHUP && renderQueue.get()) {
		renderQueue->flush();
	} else if (signal == SIGFPE) {
		SLog(EWarn, "Caught a floating-point exception!");
	}
}
#endif

class FlushThread : public Thread {
public:
	FlushThread(int timeout) : Thread("flush"), 
		m_flag(new WaitFlag()),
		m_timeout(timeout) { }

	void run() {
		while (!m_flag->get()) {
			m_flag->wait(m_timeout * 1000);
			renderQueue->flush();
		}
	}

	void quit() {
		m_flag->set(true);
		join();
	}
private:
	ref<WaitFlag> m_flag;
	int m_timeout;
};

int ubi_main(int argc, char **argv) {
	char optchar, *end_ptr = NULL;

	try {
		/* Default settings */
		int nprocs = getProcessorCount(), numParallelScenes = 1;
		std::string nodeName = getHostName(),
					networkHosts = "", destFile="";
		bool quietMode = false, progressBars = true, skipExisting = false;
		ELogLevel logLevel = EInfo;
		ref<FileResolver> fileResolver = Thread::getThread()->getFileResolver();
		bool testCaseMode = false, treatWarningsAsErrors = false;
		std::map<std::string, std::string> parameters;
		int blockSize = 32;
		int flushTimer = -1;

		if (argc < 2) {
			help();
			return 0;
		}

		optind = 1;
		/* Parse command-line arguments */
		while ((optchar = getopt(argc, argv, "a:c:D:s:j:n:o:r:b:p:qhzvtwx")) != -1) {
			switch (optchar) {
				case 'a': {
						std::vector<std::string> paths = tokenize(optarg, ";");
						for (unsigned int i=0; i<paths.size(); ++i) 
							fileResolver->addPath(paths[i]);
					}
					break;
				case 'c':
					networkHosts = networkHosts + std::string(";") + std::string(optarg);
					break;
				case 'w':
					treatWarningsAsErrors = true;
					break;
				case 'D': {
						std::vector<std::string> param = tokenize(optarg, "=");
						if (param.size() != 2)
							SLog(EError, "Invalid parameter specification \"%s\"", optarg);
						parameters[param[0]] = param[1];
					}
					break;
				case 's': {
						std::ifstream is(optarg);
						if (is.fail())
							SLog(EError, "Could not open host file!");
						std::string host;
						while (is >> host) {
							if (host.length() < 1 || host.c_str()[0] == '#')
								continue;
							networkHosts = networkHosts + std::string(";") + host;
						}
					}
					break;
				case 'n':
					nodeName = optarg;
					break;
				case 'o':
					destFile = optarg;
					break;
				case 'v':
					logLevel = EDebug;
					break;
				case 't':
					testCaseMode = true;
					break;
				case 'x':
					skipExisting = true;
					break;
				case 'p':
					nprocs = strtol(optarg, &end_ptr, 10);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the processor count!");
					break;
				case 'j':
					numParallelScenes = strtol(optarg, &end_ptr, 10);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the parallel scene count!");
					break;
				case 'r':
					flushTimer = strtol(optarg, &end_ptr, 10);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the '-r' parameter argument!");
					break;
				case 'b':
					blockSize = strtol(optarg, &end_ptr, 10);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the block size!");
					if (blockSize < 2 || blockSize > 64)
						SLog(EError, "Invalid block size (should be in the range 2-64)");
					break;
				case 'z':
					progressBars = false;
					break;
				case 'q':
					quietMode = true;
					break;
				case 'h':
				default:
					help();
					return 0;
			}
		}

		ProgressReporter::setEnabled(progressBars);

		/* Initialize OpenMP */
		Thread::initializeOpenMP(nprocs);

		/* Configure the logging subsystem */
		ref<Logger> log = Thread::getThread()->getLogger();
		log->setLogLevel(logLevel);
		log->setErrorLevel(treatWarningsAsErrors ? EWarn : EError);

		/* Disable the default appenders */
		for (size_t i=0; i<log->getAppenderCount(); ++i) {
			Appender *appender = log->getAppender(i);
			if (appender->getClass()->derivesFrom(MTS_CLASS(StreamAppender)))
				log->removeAppender(appender);
		}

		log->addAppender(new StreamAppender(formatString("mitsuba.%s.log", nodeName.c_str())));
		if (!quietMode)
			log->addAppender(new StreamAppender(&std::cout));

		SLog(EInfo, "Mitsuba version " MTS_VERSION ", Copyright (c) " MTS_YEAR " Wenzel Jakob");

		/* Configure the scheduling subsystem */
		Scheduler *scheduler = Scheduler::getInstance();
		for (int i=0; i<nprocs; ++i)
			scheduler->registerWorker(new LocalWorker(formatString("wrk%i", i)));
		std::vector<std::string> hosts = tokenize(networkHosts, ";");

		/* Establish network connections to nested servers */ 
		for (size_t i=0; i<hosts.size(); ++i) {
			const std::string &hostName = hosts[i];
			ref<Stream> stream;

			if (hostName.find("@") == std::string::npos) {
				int port = MTS_DEFAULT_PORT;
				std::vector<std::string> tokens = tokenize(hostName, ":");
				if (tokens.size() == 0 || tokens.size() > 2) {
					SLog(EError, "Invalid host specification '%s'!", hostName.c_str());
				} else if (tokens.size() == 2) {
					port = strtol(tokens[1].c_str(), &end_ptr, 10);
					if (*end_ptr != '\0')
						SLog(EError, "Invalid host specification '%s'!", hostName.c_str());
				}
				stream = new SocketStream(tokens[0], port);
			} else {
				std::string path = "~/mitsuba"; // default path if not specified
				std::vector<std::string> tokens = tokenize(hostName, "@:");
				if (tokens.size() < 2 || tokens.size() > 3) {
					SLog(EError, "Invalid host specification '%s'!", hostName.c_str());
				} else if (tokens.size() == 3) {
					path = tokens[2];
				}
				std::vector<std::string> cmdLine;
				cmdLine.push_back(formatString("bash -c 'cd %s; . setpath.sh; mtssrv -ls'", path.c_str()));
				stream = new SSHStream(tokens[0], tokens[1], cmdLine);
			}
			try {
				scheduler->registerWorker(new RemoteWorker(formatString("net%i", i), stream));
			} catch (std::runtime_error &e) {
				if (hostName.find("@") != std::string::npos) {
#if defined(WIN32)
					SLog(EWarn, "Please ensure that passwordless authentication "
						"using plink.exe and pageant.exe is enabled (see the documentation for more information)");
#else
					SLog(EWarn, "Please ensure that passwordless authentication "
						"is enabled (e.g. using ssh-agent - see the documentation for more information)");
#endif
				}
				throw e;
			}
		}

		scheduler->start();

#if !defined(WIN32)
			/* Initialize signal handlers */
			struct sigaction sa;
			sa.sa_handler = signalHandler;
			sigemptyset(&sa.sa_mask);
			sa.sa_flags = 0;
			if (sigaction(SIGHUP, &sa, NULL)) 
				SLog(EError, "Could not install a custom signal handler!");
			if (sigaction(SIGFPE, &sa, NULL)) 
				SLog(EError, "Could not install a custom signal handler!");
#endif

		/* Prepare for parsing scene descriptions */
		SAXParser* parser = new SAXParser();
		fs::path schemaPath = fileResolver->resolveAbsolute("schema/scene.xsd");

		/* Check against the 'scene.xsd' XML Schema */
		parser->setDoSchema(true);
		parser->setValidationSchemaFullChecking(true);
		parser->setValidationScheme(SAXParser::Val_Always);
		parser->setExternalNoNamespaceSchemaLocation(schemaPath.file_string().c_str());
		parser->setCalculateSrcOfs(true);

		/* Set the handler */
		SceneHandler *handler = new SceneHandler(parser, parameters);
		parser->setDoNamespaces(true);
		parser->setDocumentHandler(handler);
		parser->setErrorHandler(handler);

		renderQueue = new RenderQueue();
	
		ref<TestSupervisor> testSupervisor;

		if (testCaseMode)
			testSupervisor = new TestSupervisor(argc-optind);

		ref<FlushThread> flushThread;
		if (flushTimer > 0) {
			flushThread = new FlushThread(flushTimer);
			flushThread->start();
		}

		int jobIdx = 0;
		for (int i=optind; i<argc; ++i) {
			fs::path 
				filename = fileResolver->resolve(argv[i]),
				filePath = fs::complete(filename).parent_path(),
				baseName = fs::basename(filename);
			ref<FileResolver> frClone = fileResolver->clone();
			frClone->addPath(filePath);
			Thread::getThread()->setFileResolver(frClone);

			SLog(EInfo, "Parsing scene description from \"%s\" ..", argv[i]);

			parser->parse(filename.file_string().c_str());
			ref<Scene> scene = handler->getScene();

			scene->setSourceFile(filename);
			scene->setDestinationFile(destFile.length() > 0 ? 
				fs::path(destFile) : (filePath / baseName));
			scene->setBlockSize(blockSize);

			if (scene->destinationExists() && skipExisting)
				continue;

			ref<RenderJob> thr = new RenderJob(formatString("ren%i", jobIdx++), 
				scene, renderQueue, testSupervisor, -1, -1, -1, true,
				flushTimer > 0);
			thr->start();

			renderQueue->waitLeft(numParallelScenes-1);
		}

		/* Wait for all render processes to finish */
		renderQueue->waitLeft(0);
		if (flushThread)
			flushThread->quit();
		renderQueue = NULL;

		delete handler;
		delete parser;

		Statistics::getInstance()->printStats();

		if (testCaseMode) 
			testSupervisor->printSummary();
	} catch (const std::exception &e) {
		std::cerr << "Caught a critical exeption: " << e.what() << std::endl;
		return -1;
	} catch (...) {
		std::cerr << "Caught a critical exeption of unknown type!" << endl;
		return -1;
	}

	return 0;
}

int main(int argc, char **argv) {
	/* Initialize the core framework */
	Class::staticInitialization();
	PluginManager::staticInitialization();
	Statistics::staticInitialization();
	Thread::staticInitialization();
	Logger::staticInitialization();
	Spectrum::staticInitialization();
	Scheduler::staticInitialization();
	SHVector::staticInitialization();

#ifdef WIN32
	/* Initialize WINSOCK2 */
	WSADATA wsaData;
	if (WSAStartup(MAKEWORD(2,2), &wsaData)) 
		SLog(EError, "Could not initialize WinSock2!");
	if (LOBYTE(wsaData.wVersion) != 2 || HIBYTE(wsaData.wVersion) != 2)
		SLog(EError, "Could not find the required version of winsock.dll!");
#endif

#if !defined(WIN32)
	/* Correct number parsing on some locales (e.g. ru_RU) */
	setlocale(LC_NUMERIC, "C");
#endif

	/* Initialize Xerces-C */
	try {
		XMLPlatformUtils::Initialize();
	} catch(const XMLException &toCatch) {
		SLog(EError, "Error during Xerces initialization: %s",
			XMLString::transcode(toCatch.getMessage()));
		return -1;
	}
	
	int retval = ubi_main(argc, argv);

	XMLPlatformUtils::Terminate();

	/* Shutdown the core framework */
	SHVector::staticShutdown();
	Scheduler::staticShutdown();
	Spectrum::staticShutdown();
	Logger::staticShutdown();
	Thread::staticShutdown();
	Statistics::staticShutdown();
	PluginManager::staticShutdown();
	Class::staticShutdown();
	
#ifdef WIN32
	/* Shut down WINSOCK2 */
	WSACleanup();
#endif

	return retval;
}
