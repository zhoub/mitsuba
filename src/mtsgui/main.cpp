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

#include <QtGui/QtGui>
#include <mitsuba/core/shvector.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/appender.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/scenehandler.h>

#if defined(__OSX__)
#include <ApplicationServices/ApplicationServices.h>
#endif

#include "mainwindow.h"

#if defined(__LINUX__)
#include <X11/Xlib.h>
#endif

#if !defined(__WINDOWS__)
#include <signal.h>
#include <sys/wait.h>
#include <errno.h>
#endif

#if defined(MTS_HAS_BREAKPAD)
#if defined(__LINUX__)
#include <client/linux/handler/exception_handler.h>
#elif defined(__WINDOWS__)
#include <client/windows/handler/exception_handler.h>
#elif defined(__OSX__)
#include <client/mac/handler/exception_handler.h>
#endif
#endif

XERCES_CPP_NAMESPACE_USE

using namespace mitsuba;

MainWindow *mainWindow = NULL;

class MitsubaApplication : public QApplication {
public:
	MitsubaApplication(int &argc, char **argv) : QApplication(argc, argv) {
	}

	bool event(QEvent *event) {
		switch (event->type()) {
#if defined(__OSX__)
			case QEvent::Quit:
				quit();
				return true;
#endif
			case QEvent::FileOpen:
				if (mainWindow != NULL)
					mainWindow->loadFile(static_cast<QFileOpenEvent *>(event)->file());
				return true;
			default:
				return QApplication::event(event);
		}
	}

	bool notify(QObject *rec, QEvent *e) {
		try {
			return QApplication::notify(rec, e);
		} catch (const std::exception &e) {
			SLog(EWarn, "Caught exception: %s", e.what());
			QMessageBox::critical(NULL, tr("Critical exception"),
				e.what(), QMessageBox::Ok);
			return false;
		}
	}
};

/* Collect zombie processes */
#if !defined(__WINDOWS__)
void collect_zombies(int s) {
	while (waitpid(-1, NULL, WNOHANG) > 0);
}
#endif

#if defined(MTS_HAS_BREAKPAD) && defined(__LINUX__)
static bool dumpCallbackLinux(
		const char* dump_path,
		const char* minidump_id,
		void* context, bool succeeded) {
	SLog(EInfo, "A critical application error occurred! A mini-dump file "
			"has been created in %s/%s.dmp. Please submit this file "
			" and the responsible scene to the Mitsuba bug tracker at "
			"https://www.mitsuba-renderer.org/tracker.", dump_path, minidump_id);
	return succeeded;
}
#endif

int main(int argc, char *argv[]) {
	int retval;

#if defined(MTS_HAS_BREAKPAD)
	/* Set up Google Breakpad, a crash handler that can be 
	   used to get useful backtraces from end-users. */
	
	#if defined(__LINUX__)
		google_breakpad::ExceptionHandler eh(
			"/tmp", NULL, dumpCallbackLinux, NULL, true);
	#endif
#endif

	/* Initialize the core framework */
	Class::staticInitialization();
	PluginManager::staticInitialization();
	Statistics::staticInitialization();
	Thread::staticInitialization();
	Thread::initializeOpenMP(getProcessorCount());
	Logger::staticInitialization();
	Spectrum::staticInitialization();
	Scheduler::staticInitialization();
	SHVector::staticInitialization();
	SceneHandler::staticInitialization();

#if defined(__LINUX__)
	XInitThreads();
#endif

#if defined(__OSX__)
	MTS_AUTORELEASE_BEGIN()
	/* Required for the mouse relocation in GLWidget */
	CGSetLocalEventsSuppressionInterval(0.0f);
	MTS_AUTORELEASE_END() 
#endif

#if defined(__WINDOWS)
	/* Initialize WINSOCK2 */
	WSADATA wsaData;
	if (WSAStartup(MAKEWORD(2,2), &wsaData)) 
		SLog(EError, "Could not initialize WinSock2!");
	if (LOBYTE(wsaData.wVersion) != 2 || HIBYTE(wsaData.wVersion) != 2)
		SLog(EError, "Could not find the required version of winsock.dll!");
#endif

#if !defined(__WINDOWS__)
	/* Avoid zombies processes when running the server */
	struct sigaction sa;
	sa.sa_handler = collect_zombies;
	sigemptyset(&sa.sa_mask);
	sa.sa_flags = SA_RESTART;

	if (sigaction(SIGCHLD, &sa, NULL) == -1)
		SLog(EWarn, "Error in sigaction(): %s!", strerror(errno));
#endif

	qRegisterMetaType<ELogLevel>("ELogLevel");
	qRegisterMetaType<fs::path>("fs::path");

	MitsubaApplication app(argc, argv);
	try {
		QFile stylesheet(":/resources/stylesheet.css");

		if (!stylesheet.open(QFile::ReadOnly)) {
			QMessageBox::critical(NULL, "Internal error", "Could not open stylesheet!");
			exit(-1);
		}
		app.setStyleSheet(QTextStream(&stylesheet).readAll().toAscii());

#if defined(__OSX__)
		app.setAttribute(Qt::AA_DontShowIconsInMenus); 
#endif
		/* Disable the default appenders */
		ref<Logger> logger = Thread::getThread()->getLogger();
		for (size_t i=0; i<logger->getAppenderCount(); ++i) {
			Appender *appender = logger->getAppender(i);
			if (appender->getClass()->derivesFrom(MTS_CLASS(StreamAppender)))
				logger->removeAppender(appender);
		}

#if defined(__OSX__)
		/* Create a log file inside the application bundle */
		MTS_AUTORELEASE_BEGIN() 
		logger->addAppender(new StreamAppender(formatString("%s/mitsuba.%s.log", 
			__mts_bundlepath().c_str(), getHostName().c_str())));
		MTS_AUTORELEASE_END() 
#else
		/* Create a log file inside the current working directory */
		logger->addAppender(new StreamAppender(formatString("mitsuba.%s.log", getHostName().c_str())));
#endif

#if !defined(__WINDOWS__)
		/* Correct number parsing on some locales (e.g. ru_RU) */
		setlocale(LC_NUMERIC, "C");
#endif

		mainWindow = new MainWindow();
		mainWindow->initWorkers();
		retval = app.exec();
		delete mainWindow;
	} catch (const std::exception &e) {
		SLog(EWarn, "Critical exception during startup: %s", e.what());
		QMessageBox::critical(NULL, QString("Critical exception"),
			e.what(), QMessageBox::Ok);
		retval = -1;
	}
	Statistics::getInstance()->printStats();


#if defined(__WINDOWS__)
	/* Shut down WINSOCK2 */
	WSACleanup();
#endif

	/* Shutdown the core framework */
	SceneHandler::staticShutdown();
	SHVector::staticShutdown();
	Scheduler::staticShutdown();
	Spectrum::staticShutdown();
	Logger::staticShutdown();
	Thread::staticShutdown();
	Statistics::staticShutdown();
	PluginManager::staticShutdown();
	Class::staticShutdown();

	return retval;
}
