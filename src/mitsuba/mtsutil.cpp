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

#include <mitsuba/core/platform.h>
#include <xercesc/parsers/SAXParser.hpp>
#include <mitsuba/core/sched_remote.h>
#include <mitsuba/core/sstream.h>
#include <mitsuba/core/sshstream.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/appender.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/renderjob.h>
#include <mitsuba/render/scenehandler.h>
#include <boost/algorithm/string.hpp>
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
	cout <<  "Usage: mtsutil [mtsutil options] <utility name> [arguments]" << endl;
	cout <<  "Options/Arguments:" << endl;
	cout <<  "   -h          Display this help text" << endl << endl;
	cout <<  "   -a p1;p2;.. Add one or more entries to the resource search path" << endl << endl;
	cout <<  "   -p count    Override the detected number of processors. Useful for reducing" << endl;
	cout <<  "               the load or creating scheduling-only nodes in conjunction with"  << endl;
	cout <<  "               the -c and -s parameters, e.g. -p 0 -c host1;host2;host3,..." << endl << endl;
	cout <<  "   -q          Quiet mode - do not print any log messages to stdout" << endl << endl;
	cout <<  "   -c hosts    Network processing: connect to mtssrv instances over a network." << endl;
	cout <<  "               Requires a semicolon-separated list of host names of the form" << endl;
	cout <<  "                       host.domain[:port] for a direct connection" << endl;
	cout <<  "                 or" << endl;
	cout <<  "                       user@host.domain[:path] for a SSH connection (where" << endl;
	cout <<  "                       \"path\" denotes the place where Mitsuba is checked" << endl;
	cout <<  "                       out -- by default, \"~/mitsuba\" is used)" << endl << endl;
	cout <<  "   -s file     Connect to additional Mitsuba servers specified in a file" << endl;
	cout <<  "               with one name per line (same format as in -c)" << endl<< endl;
	cout <<  "   -n name     Assign a node name to this instance (Default: host name)" << endl << endl;
	cout <<  "   -t          Execute all testcases" << endl << endl;
	cout <<  "   -v          Be more verbose" << endl << endl;

	FileResolver *fileResolver = Thread::getThread()->getFileResolver();
	std::ostringstream utilities, testcases;

	testcases << "The following testcases are available:" << endl << endl;
	utilities << endl << "The following utilities are available:" << endl << endl;

	std::vector<fs::path> dirPaths = fileResolver->resolveAll("plugins");
	std::set<std::string> seen;

	for (size_t i=0; i<dirPaths.size(); ++i) {
		fs::path dirPath = fs::complete(dirPaths[i]);

		if (!fs::exists(dirPath) || !fs::is_directory(dirPath))
			break;

		fs::directory_iterator end, it(dirPath);

		for (; it != end; ++it) {
			if (!fs::is_regular_file(it->status()))
				continue;
			std::string extension(boost::to_lower_copy(it->path().extension()));
#if defined(WIN32)
			if (extension != ".dll")
				continue;
#elif defined(__OSX__)
			if (extension != ".dylib")
				continue;
#elif defined(__LINUX__)
			if (extension != ".so")
				continue;
#else
#error Unknown operating system!
#endif
			std::string shortName = it->path().stem();
			if (seen.find(shortName) != seen.end())
				continue;
			seen.insert(shortName);
			Plugin utility(shortName, it->path());
			if (!utility.isUtility())
				continue;
			if (boost::starts_with(shortName, "test_")) {
				testcases << "\t" << shortName;
				for (int i=0; i<22-(int) shortName.length(); ++i)
					testcases << ' ';
				testcases << utility.getDescription() << endl;
			} else {
				utilities << "\t" << shortName;
				for (int i=0; i<22-(int) shortName.length(); ++i)
					utilities << ' ';
				utilities << utility.getDescription() << endl;
			}
		}
	}

	cout << testcases.str() << utilities.str();
}


int mtsutil(int argc, char **argv) {
	char optchar, *end_ptr = NULL;

	try {
		/* Default settings */
		int nprocs = getProcessorCount();
		std::string nodeName = getHostName(),
					networkHosts = "", destFile="";
		bool quietMode = false;
		ELogLevel logLevel = EInfo;
		FileResolver *fileResolver = Thread::getThread()->getFileResolver();
		bool testCaseMode = false;

		if (argc < 2) {
			help();
			return 0;
		}

		optind = 1;
		/* Parse command-line arguments */
		while ((optchar = getopt(argc, argv, "+a:c:s:n:p:qhvt")) != -1) {
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
				case 't':
					testCaseMode = true;
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
				case 'v':
					logLevel = EDebug;
					break;
				case 'p':
					nprocs = strtol(optarg, &end_ptr, 10);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the processor count!");
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

		/* Configure the logging subsystem */
		ref<Logger> log = Thread::getThread()->getLogger();
		log->setLogLevel(logLevel);
	
		/* Initialize OpenMP */
		Thread::initializeOpenMP(nprocs);
		
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

		if (testCaseMode) {
			std::vector<fs::path> dirPaths = fileResolver->resolveAll("plugins");
			std::set<std::string> seen;
			int executed = 0, succeeded = 0;
		
			for (size_t i=0; i<dirPaths.size(); ++i) {
				fs::path dirPath = fs::complete(dirPaths[i]);

				if (!fs::exists(dirPath) || !fs::is_directory(dirPath))
					break;

				fs::directory_iterator end, it(dirPath);

				for (; it != end; ++it) {
					if (!fs::is_regular_file(it->status()))
						continue;
					std::string extension(boost::to_lower_copy(it->path().extension()));
#if defined(WIN32)
					if (extension != ".dll")
						continue;
#elif defined(__OSX__)
					if (extension != ".dylib")
						continue;
#elif defined(__LINUX__)
					if (extension != ".so")
						continue;
#else
#error Unknown operating system!
#endif
					std::string shortName = it->path().stem();
					if (seen.find(shortName) != seen.end() || !boost::starts_with(shortName, "test_"))
						continue;
					seen.insert(shortName);
					Plugin plugin(shortName, it->path());
					if (!plugin.isUtility())
						continue;

					ref<Utility> utility = plugin.createUtility();

					TestCase *testCase = static_cast<TestCase *>(utility.get());
					if (!utility->getClass()->derivesFrom(MTS_CLASS(TestCase)))
						SLog(EError, "This is not a test case!");

					if (testCase->run(argc-optind, argv+optind) != 0)
						SLog(EError, "Testcase unexpectedly returned with a nonzero value.");

					executed += testCase->getExecuted();
					succeeded += testCase->getSucceeded();
				}
			}

			SLog(EInfo, "Ran %i tests, %i succeeded, %i failed.", executed, succeeded, executed-succeeded);
		} else {
			if (argc <= optind) {
				std::cerr << "A utility name must be supplied!" << endl;
				return -1;
			}
			fs::path pluginName(argv[optind]);

			/* Build the full plugin file name */
#if defined(WIN32)
			pluginName.replace_extension(".dll");
#elif defined(__OSX__)
			pluginName.replace_extension(".dylib");
#elif defined(__LINUX__)
			pluginName.replace_extension(".so");
#else
#error Unknown operating system!
#endif
			fs::path fullName = fileResolver->resolve(fs::path("plugins") / pluginName);

			if (!fs::exists(fullName)) {
				/* Plugin not found! */
				SLog(EError, "Utility \"%s\" not found (run \"mtsutil\" without arguments to "
					"see a list of available utilities)", fullName.file_string().c_str());
			}

			SLog(EInfo, "Loading utility \"%s\" ..", argv[optind]);
			Plugin *plugin = new Plugin(argv[optind], fullName);
			if (!plugin->isUtility())
				SLog(EError, "This plugin does not implement the 'Utility' interface!");
			Statistics::getInstance()->logPlugin(argv[optind], plugin->getDescription());
		
			ref<Utility> utility = plugin->createUtility();

			int retval = utility->run(argc-optind, argv+optind);
			scheduler->pause();
			utility = NULL;
			scheduler->stop();
			delete plugin;
			return retval;
		}
	} catch (const std::exception &e) {
		std::cerr << "Caught a critical exeption: " << e.what() << std::endl;
	} catch (...) {
		std::cerr << "Caught a critical exeption of unknown type!" << endl;
	}

	return 0;
}

int ubi_main(int argc, char **argv) {
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
	
	int retval = mtsutil(argc, argv);

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

#if !defined(__OSX__)
int main(int argc, char **argv) {
	return ubi_main(argc, argv);
}
#endif


