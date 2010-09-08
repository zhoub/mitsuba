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

#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "rendersettingsdlg.h"
#include "previewsettingsdlg.h"
#include "programsettingsdlg.h"
#include "sceneloader.h"
#include "logwidget.h"
#include "navdlg.h"
#include "aboutdlg.h"
#include "importdlg.h"
#include "server.h"
#include "save.h"
#include <QtNetwork>
#include <mitsuba/core/sched_remote.h>
#include <mitsuba/core/sstream.h>
#include <mitsuba/core/sshstream.h>
#include <mitsuba/core/plugin.h>

#if !defined(WIN32)
#include <QX11Info>
#include <pwd.h>
#endif

#if defined(__OSX__)
#include "previewsettingsdlg_cocoa.h"
#else
#include "previewsettingsdlg.h"
#endif

class ProgramVersion {
public:
	ProgramVersion(const QString &versionString) {
		QStringList sl;
		sl = versionString.trimmed().split('.');
		SAssert(sl.size() == 3);
		major = sl[0].toInt();
		minor = sl[1].toInt();
		release = sl[2].toInt();
	}

	inline bool operator<(const ProgramVersion &other) const {
		if (major < other.major)
			return true;
		if (major > other.major)
			return false;
		if (minor < other.minor)
			return true;
		if (minor > other.minor)
			return false;
		if (release < other.release)
			return true;
		return false;
	}

	QString toString() const {
		return QString("%1.%2.%3").arg(major).arg(minor).arg(release);
	}
private:
	int major;
	int minor;
	int release;
};

static int localWorkerCtr = 0, remoteWorkerCtr = 0;

MainWindow::MainWindow(QWidget *parent) :
	QMainWindow(parent), ui(new Ui::MainWindow), 
	m_networkReply(NULL), m_activeWindowHack(false) {
	Logger *logger = Thread::getThread()->getLogger();

#if defined(__OSX__)
	m_previewSettings = NULL;
#endif

	QSettings settings("mitsuba-renderer.org", "qtgui");
	logger->setLogLevel((ELogLevel) settings.value("verbosity", EDebug).toInt());

	m_logWidget = new LogWidget(NULL);

	m_consoleAppender = new QConsoleAppender();
	logger->addAppender(m_consoleAppender);
	connect(m_consoleAppender, SIGNAL(textMessage(ELogLevel, const QString &)), 
		m_logWidget, SLOT(onTextMessage(ELogLevel, const QString &)), Qt::QueuedConnection);

	SLog(EInfo, "Mitsuba version " MTS_VERSION ", Copyright (c) " MTS_YEAR " Wenzel Jakob");

	m_currentChild = NULL;
	ui->setupUi(this);
	isActive();
	m_lastTab = NULL;
	ui->glView->setScrollBars(ui->hScrollBar, ui->vScrollBar);
	ui->hScrollBar->setVisible(false);
	ui->vScrollBar->setVisible(false);
	ui->actionUpdateCheck->setMenuRole(QAction::ApplicationSpecificRole);
	m_progressWidget = new QWidget(centralWidget());
	m_progressLabel = new QLabel(m_progressWidget);
	m_progress = new QProgressBar(m_progressWidget);
	m_progressWidget->setObjectName("progressWidget");
	QHBoxLayout *hlayout = new QHBoxLayout();
	hlayout->addWidget(m_progressLabel);
	hlayout->addWidget(m_progress);
	m_progressWidget->setLayout(hlayout);
#if defined(__OSX__)
	m_progressWidget->setStyleSheet("QWidget#progressWidget {background: qlineargradient(x1:0, y1:0, x2:0, y2:1, "
		"stop:0 palette(dark), stop: 1 palette(mid)); border-top: 1px solid palette(mid); margin: 0px; spacing: 15px; }");
	m_progress->setMaximumSize(SHRT_MAX,18);
	m_progress->setMinimumSize(10,18);
	hlayout->setContentsMargins(10, 8, 10, 3);
	m_progress->setAttribute(Qt::WA_MacSmallSize, true);
	m_progressLabel->setAttribute(Qt::WA_MacSmallSize, true);
#endif
	m_serverWidget = NULL;
	m_contextIndex = -1;

	for (int i = 0; i < MAX_RECENT_FILES; ++i) {
		m_actRecent[i] = new QAction(this);
		m_actRecent[i]->setVisible(false);
		connect(m_actRecent[i], SIGNAL(triggered()),
			this, SLOT(onOpenRecent()));
		ui->menuOpen_Recent->addAction(m_actRecent[i]);
	}

	m_clearRecent = new QAction(tr("Clear Menu"), this);
	connect(m_clearRecent, SIGNAL(triggered()),
		this, SLOT(onClearRecent()));
	ui->menuOpen_Recent->addAction(m_clearRecent);

#if defined(__OSX__)
	/* Make this the default menu bar */
	ui->menuBar->setParent(NULL);
#endif

	updateRecentFileActions();

	ui->tabBar->setDocumentMode(true);
	ui->tabBar->setTabsClosable(true);
	ui->tabBar->setMovable(true);
	ui->tabBar->setContextMenuPolicy(Qt::CustomContextMenu);
	on_tabBar_currentChanged(-1);

	connect(ui->glView, SIGNAL(quit()), this, SLOT(on_actionExit_triggered()));
	connect(ui->glView, SIGNAL(beginRendering()), this, SLOT(on_actionRender_triggered()));
	connect(ui->glView, SIGNAL(stopRendering()), this, SLOT(updateUI()));
	connect(ui->glView, SIGNAL(statusMessage(const QString &)), this, SLOT(onStatusMessage(const QString &)));

	/* Load defaults from app settings file */
	ui->glView->setInvertMouse(settings.value("invertMouse", false).toBool());
	ui->glView->setMouseSensitivity(settings.value("mouseSensitivity", 3).toInt());
	ui->glView->setNavigationMode((ENavigationMode) settings.value("navigationMode", 
		EFlythroughFixedYaw).toInt());
	m_searchPaths = settings.value("searchPaths", QStringList()).toStringList();
	m_blockSize = settings.value("blockSize", 32).toInt();
	m_listenPort = settings.value("listenPort", MTS_DEFAULT_PORT).toInt();
	m_nodeName = settings.value("nodeName", getFQDN().c_str()).toString();

	m_renderQueue = new RenderQueue();
	m_renderListener = new QRenderListener();
	m_renderQueue->registerListener(m_renderListener);

	connect(m_renderListener, SIGNAL(jobFinished(const RenderJob *, bool)), 
		this, SLOT(onJobFinished(const RenderJob *, bool)), Qt::QueuedConnection);
	connect(m_renderListener, SIGNAL(refresh(const RenderJob *)), 
		this, SLOT(onRefresh(const RenderJob *)), Qt::QueuedConnection);
	connect(m_renderListener, SIGNAL(workEnd(const RenderJob *, const ImageBlock *)), 
		this, SLOT(onWorkEnd(const RenderJob *, const ImageBlock *)), Qt::DirectConnection);
	connect(m_renderListener, SIGNAL(workBegin(const RenderJob *, const RectangularWorkUnit *, int)),
        this, SLOT(onWorkBegin(const RenderJob *, const RectangularWorkUnit *, int)), Qt::DirectConnection);
	connect(m_consoleAppender, 
		SIGNAL(progressMessage(const RenderJob *, const QString &, float, const QString &)), 
		this, SLOT(onProgressMessage(const RenderJob *, const QString &, float, const QString &)), 
		Qt::QueuedConnection);
	connect(this, SIGNAL(updateView()), ui->glView, SLOT(onUpdateView()));
	
	QPoint windowPos;
	if (settings.contains("pos")) {
		windowPos = settings.value("pos").toPoint();
	} else {
		QDesktopWidget *desktop = QApplication::desktop();
		windowPos = QPoint((desktop->width() - width()) / 2, (desktop->height() - height())/2);
	}

#if defined(__OSX__)
	ui->toolBar->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
	QToolButton *previewButton = static_cast<QToolButton *>(ui->toolBar->widgetForAction(ui->actionPreviewSettings));
	previewButton->setStyleSheet("margin-left: -5px; margin-right:-5px");
	
	/* Weird Qt/OSX bug -- moving while a window while it is invisible causes 
	   it to appear move up by 65 pixels (this is related to the unified toolbar) */
	move(windowPos + QPoint(0, 65));
#else
	move(windowPos);
#endif
	show();
	/* Move again just to be sure.. */
	move(windowPos);

	updateUI();
	adjustSize();

	m_networkManager = new QNetworkAccessManager(this);

#if defined(__OSX__)
	/* Submit crash reports on OSX */
	QDir crashDir = QDir::home();
	crashDir.cd("Library/Logs/CrashReporter");
	QFileInfoList crashReports = crashDir.entryInfoList(QStringList("qtgui_*"), 
			QDir::Files, QDir::Name);

	if (crashReports.size() > 0) {
		if (QMessageBox::question(this, tr("Crash reports detected"),
			QString("<p>It appears that Mitsuba has crashed on a previous run!</p><p>"
				"%1 crash reports were found. Would you like to submit them?").arg(crashReports.size()),
				QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes) {
			bool failure = false;
			for (int i=0; i<crashReports.size(); ++i) {
				QFile file(crashReports[i].absoluteFilePath());
				if (!file.open(QIODevice::ReadOnly)) {
					failure = true;
					QMessageBox::critical(this, tr("Unable to read crash report"),
						tr("Unable to read crash report -- please check the file permissions in ~/Library/Logs/CrashReporter."), QMessageBox::Ok);
					break;
				}

				QString username("unknown");
				uid_t uid = getuid();
				struct passwd *info = getpwuid(uid);
				if (info) 
					username = QString(info->pw_name);

				QString boundary("5da85133908e2"); // some arbitrary value
				QByteArray data(QString("--" + boundary + "\r\n").toAscii());
				data.append(QString("Content-Disposition: form-data; name=\"bugreport\"; filename=\"%1\"\r\n").arg(file.fileName()));
				data.append("Content-Type: application/octet-stream\r\n\r\n");
				QString header = QString("Bug report from machine \"%1\", user \"%2\", filename \"%3\""
					", Mitsuba version " MTS_VERSION).arg(getFQDN().c_str()).arg(username).arg(file.fileName());
				data.append(header + "\r\n");
				for (int j=0; j<header.length(); j++)
					data.append('=');
				data.append("\r\n\r\n");
				data.append(file.readAll());
				data.append("\r\n--" + boundary + "--\r\n");

				QNetworkRequest reqHeader;
				reqHeader.setHeader(QNetworkRequest::ContentTypeHeader, QString("multipart/form-data; boundary=\"%1\"").arg(boundary));
				reqHeader.setHeader(QNetworkRequest::ContentLengthHeader, QString::number(data.length()));
				reqHeader.setUrl(QUrl("http://www.mitsuba-renderer.org/bugreporter.php"));
				m_bugStatus = 0;
				QNetworkReply *reply = m_networkManager->post(reqHeader, data);
				connect(reply, SIGNAL(error(QNetworkReply::NetworkError)), this, SLOT(onBugReportError()));
				connect(reply, SIGNAL(finished()), this, SLOT(onBugReportSubmitted()));
				while (m_bugStatus == 0) {
					qApp->processEvents();
					Thread::sleep(50);
				}
				if (m_bugStatus != 1) {
					failure = true;
					QMessageBox::critical(this, tr("Unable to submit"),
						tr("Unable to submit crash report -- are you connected to the internet?"), QMessageBox::Ok);
					break;
				}
				reply->deleteLater();
				if (!file.remove()) {
					failure = true;
					QMessageBox::critical(this, tr("Unable to delete submitted crash report"),
						tr("Unable to delete a submitted crash report -- please check the file permissions in ~/Library/Logs/CrashReporter."), QMessageBox::Ok);
					break;
				}
			}
			if (!failure)
				QMessageBox::information(this, tr("Crash reporter"),
					tr("All crash reports have been submitted. Thank you!"),
					QMessageBox::Ok);
		}
	}
#endif

	if (ui->glView->isUsingSoftwareFallback())
		QMessageBox::warning(this, tr("Insufficient OpenGL capabilities"),
			ui->glView->getErrorString(), QMessageBox::Ok);

	connect(m_networkManager, SIGNAL(finished(QNetworkReply *)), this, SLOT(onNetworkFinished(QNetworkReply *)));
	m_checkForUpdates = settings.value("checkForUpdates", true).toBool();
	if (m_checkForUpdates)
		checkForUpdates(false);
}

MainWindow::~MainWindow() {
	m_renderQueue->unregisterListener(m_renderListener);
	ref<Scheduler> scheduler = Scheduler::getInstance();
	scheduler->pause();
	for (int i=0; i<m_connections.size(); ++i) {
		ServerConnection &c = m_connections[i];
		scheduler->unregisterWorker(c.worker);
	}
	if (m_networkReply)
		m_networkReply->abort();
#if defined(__OSX__)
	delete ui->menuBar;
#endif
    delete ui;
}

void MainWindow::initWorkers() {
	QSettings settings("mitsuba-renderer.org", "qtgui");
	ref<Scheduler> scheduler = Scheduler::getInstance();
	int localWorkerCount = settings.value("localWorkers", getProcessorCount()).toInt();
	for (int i=0; i<localWorkerCount; ++i)
		scheduler->registerWorker(new LocalWorker(formatString("wrk%i", localWorkerCtr++)));

	int networkConnections = 0;
	QList<QVariant> connectionData = settings.value("connections").toList();
	if (connectionData.size() > 0) {
		QDialog *dialog = new NonClosableDialog(this);
		dialog->setWindowModality(Qt::WindowModal);
		QVBoxLayout *layout = new QVBoxLayout(dialog);
		QLabel *label = new QLabel(tr("Establishing network connections .."), dialog);
		label->setAlignment(Qt::AlignCenter);
		layout->addWidget(label);
		QProgressBar *progressBar = new QProgressBar(dialog);
		progressBar->setTextVisible(false);
		dialog->resize(200, 50);
		layout->addWidget(progressBar);
		progressBar->setTextVisible(false);
		// weird, Qt/Win needs this to get a busy indicator
		progressBar->setValue(1);
		progressBar->setRange(0, 0);
		dialog->show();

		for (int i=0; i<connectionData.size(); ++i) {
			ServerConnection conn;
			conn.fromByteArray(connectionData[i].toByteArray());
			for (int j=0; j<10; ++j)
				qApp->processEvents();
			if (conn.createWorker(this)) {
				++networkConnections;
				conn.isRegistered = true;
				scheduler->registerWorker(conn.worker);
				m_connections.append(conn);
			}
		}

		/* Update, this removes connection failures */
		connectionData.clear();
		for (int i=0; i<m_connections.size(); ++i)
			connectionData.append(m_connections[i].toByteArray());
		settings.setValue("connections", connectionData);
	
		dialog->hide();
		delete dialog;
	}

	if (networkConnections + localWorkerCount == 0) {
		QMessageBox::warning(this, tr("Scheduler warning"),
			tr("There must be at least one worker thread -- forcing creation of one."),
			QMessageBox::Ok);
		scheduler->registerWorker(new LocalWorker(formatString("wrk%i", localWorkerCtr++)));
	}
	
	QStringList args = qApp->arguments();
	for (int i=1; i<args.count(); ++i)
		loadFile(args[i]);

	scheduler->start();
	raise();
}


void MainWindow::adjustSize() {
	/* Like QWidget::adjustSize, but clamps the size to 4/5ths of the screen */
	ensurePolished();
    QSize s = sizeHint();

#if defined(Q_WS_X11)
        QRect screen = QApplication::desktop()->screenGeometry(x11Info().screen());
#else // all others
        QRect screen = QApplication::desktop()->screenGeometry(pos());
#endif

	s.setWidth(qMin(s.width(), screen.width()*4/5));
	s.setHeight(qMin(s.height(), screen.height()*4/5));
	layout()->activate();

    if (s.isValid())
        resize(s);
}

void MainWindow::checkForUpdates(bool notifyIfNone) {
	m_notifyIfNoUpdates = notifyIfNone;
	m_networkReply = m_networkManager->get(QNetworkRequest(QUrl("http://www.mitsuba-renderer.org/version")));
}

void MainWindow::onNetworkFinished(QNetworkReply *reply) {
	if (reply->error() == QNetworkReply::NoError) {
		ProgramVersion remote(QString(reply->readAll()));
		ProgramVersion local(MTS_VERSION);

		if (local < remote)
			QMessageBox::information(this, tr("New version available"),
				QString("<p>A new version of Mitsuba is available!</p><p>"
					"You have version <b>%1</b> and the most recent release is version <b>%2</b>. "
					"Please visit the <a href=\"http://www.mitsuba-renderer.org\">Mitsuba website</a> to download it.</p>")
					.arg(local.toString()).arg(remote.toString()),
					QMessageBox::Ok);
		else if (m_notifyIfNoUpdates)
			QMessageBox::information(this, tr("Installed version is current"),
				QString("<p>You're up to date!</p>"
					"<p>Mitsuba <b>%1</b> is still the newest version available.</p>")
					.arg(local.toString()), QMessageBox::Ok);
	} else {
		if (m_notifyIfNoUpdates)
			QMessageBox::warning(this, tr("Unable to determine the newest version"),
				QString("<p>Unable to determine the newest Mitsuba version.</p><p>"
					"Perhaps you are not connected to the Internet?</p>"),
					QMessageBox::Ok);
	}
}

void MainWindow::onBugReportError() {
	m_bugStatus = 2;
}

void MainWindow::onBugReportSubmitted() {
	if (m_bugStatus == 0)
		m_bugStatus = 1;
}

void MainWindow::on_actionImport_triggered() {
#if defined(MTS_HAS_COLLADA)
	ref<FileResolver> resolver = FileResolver::getInstance();
	ref<FileResolver> newResolver = resolver->clone();
	for (int i=0; i<m_searchPaths.size(); ++i)
		newResolver->addPath(m_searchPaths[i].toStdString());

	ImportDialog *dialog = new ImportDialog(this, newResolver);
	dialog->setAttribute(Qt::WA_DeleteOnClose);
	connect(dialog, SIGNAL(finished(int)), this, SLOT(onImportDialogClose(int)));
	m_currentChild = dialog;
	// prevent a tab drawing artifact on Qt/OSX
	m_activeWindowHack = true;
	dialog->show();
	qApp->processEvents();
	m_activeWindowHack = false;
#else
	QMessageBox::critical(this, tr("Importer disabled"),
		tr("The importer is disabled in this build. To use it, you will need "
		"to install COLLADA-DOM and recompile Mitsuba -- please see the "
		"documentation for more details."), 
		QMessageBox::Ok);
#endif
}

void MainWindow::onImportDialogClose(int reason) {
	m_currentChild = NULL;
}

void MainWindow::on_actionDuplicateTab_triggered() {
	int currentIndex = ui->tabBar->currentIndex();
	if (m_contextIndex != -1)
		currentIndex = m_contextIndex;
	SceneContext *currentContext = m_context[currentIndex];
	SceneContext *newContext = new SceneContext(currentContext);

	m_contextMutex.lock();
	m_context.append(newContext);
	m_contextMutex.unlock();
	ui->tabBar->addTab(newContext->shortName);

	int index = ui->tabBar->count()-1;
	if (ui->tabBar->currentIndex() != index)
		ui->tabBar->setCurrentIndex(index);
}

void MainWindow::on_actionUpdateCheck_triggered() {
	checkForUpdates(true);
}

void MainWindow::changeEvent(QEvent *e) {
    QMainWindow::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui->retranslateUi(this);
        break;
    default:
        break;
    }
}
	
SceneContext *MainWindow::getContext(const RenderJob *job, bool failIfNotFound) {
	m_contextMutex.lock();
	for (int i=0; i<m_context.size(); ++i) {
		if (m_context[i]->renderJob == job) {
			m_contextMutex.unlock();
			return m_context[i];
		}
	}
	m_contextMutex.unlock();
	if (failIfNotFound)
		SLog(EError, "Internal error: could not find render context!");
	return NULL;
}

void MainWindow::onProgressMessage(const RenderJob *job, const QString &name, 
	float progress, const QString &eta) {
	SceneContext *context = getContext(job, false);
	if (context == NULL)
		return;
	context->eta = eta;
	context->progress = progress;
	context->progressName = name + ": ";
	updateUI();
}

void MainWindow::on_actionOpen_triggered() {
	QFileDialog *dialog = new QFileDialog(this, Qt::Sheet);
	dialog->setNameFilter(tr("Mitsuba scenes (*.xml);;EXR images (*.exr)"));
	dialog->setAttribute(Qt::WA_DeleteOnClose);
	dialog->setAcceptMode(QFileDialog::AcceptOpen);
	dialog->setViewMode(QFileDialog::Detail);
	dialog->setWindowModality(Qt::WindowModal);
	QSettings settings("mitsuba-renderer.org", "qtgui");
	dialog->restoreState(settings.value("fileDialogState").toByteArray());
	connect(dialog, SIGNAL(finished(int)), this, SLOT(onOpenDialogClose(int)));
	m_currentChild = dialog;
	// prevent a tab drawing artifact on Qt/OSX
	m_activeWindowHack = true;
	dialog->show();
	qApp->processEvents();
	m_activeWindowHack = false;
}

void MainWindow::onOpenDialogClose(int reason) {
	QSettings settings("mitsuba-renderer.org", "qtgui");
	QFileDialog *dialog = static_cast<QFileDialog *>(sender());
	m_currentChild = NULL;
	if (reason == QDialog::Accepted) {
		QStringList fileNames = dialog->selectedFiles();
		settings.setValue("fileDialogState", dialog->saveState());
		for (int i=0; i<fileNames.size(); ++i)
			loadFile(fileNames[i]);
	}
}

void MainWindow::on_actionExit_triggered() {
	qApp->closeAllWindows();
}

void MainWindow::onOpenRecent() {
	QAction *action = qobject_cast<QAction *>(sender());
	if (action)
		loadFile(action->data().toString());
}

void MainWindow::onClearRecent() {
	QSettings settings("mitsuba-renderer.org", "qtgui");
	for (int j = 0; j < MAX_RECENT_FILES; ++j)
		m_actRecent[j]->setVisible(false);
	settings.setValue("recentFileList", QStringList());
}

SceneContext *MainWindow::loadScene(const QString &qFileName) {
	ref<FileResolver> resolver = FileResolver::getInstance();
	std::string filename = resolver->resolveAbsolute(qFileName.toStdString());
	std::string filePath = resolver->getParentDirectory(filename);
	ref<FileResolver> newResolver = resolver->clone();
	if (!newResolver->contains(filePath))
		newResolver->addPath(filePath);
	for (int i=0; i<m_searchPaths.size(); ++i)
		newResolver->addPath(m_searchPaths[i].toStdString());
	newResolver->setCurrentDirectoryFromFile(filename);

	ref<SceneLoader> loadingThread = new SceneLoader(newResolver, filename);
	loadingThread->start();

	QDialog *dialog = new NonClosableDialog(this);
	dialog->setWindowModality(Qt::WindowModal);
	dialog->setWindowTitle("Loading ..");
	QVBoxLayout *layout = new QVBoxLayout(dialog);
	QProgressBar *progressBar = new QProgressBar(dialog);
	dialog->resize(200, 50);
	layout->addWidget(progressBar);
	progressBar->setTextVisible(false);
	// weird, Qt/Win needs this to get a busy indicator
	progressBar->setValue(1);
	progressBar->setRange(0, 0);
	dialog->show();
	progressBar->show();

	while (loadingThread->isRunning()) {
		QCoreApplication::processEvents();
		loadingThread->wait(20);
	}
	loadingThread->join();

	dialog->hide();
	delete dialog;

	SceneContext *result = loadingThread->getResult();
	if (result == NULL) {
		QMessageBox::critical(this, tr("Unable to load %1").arg(qFileName),
			QString(loadingThread->getError().c_str()),
			QMessageBox::Ok);
	}

	return result;
}

void MainWindow::loadFile(QString filename) {
	QFileInfo fInfo(filename);
	fInfo.makeAbsolute();

	SceneContext *context = loadScene(filename);
	if (context == NULL) 
		return;
	m_contextMutex.lock();
	m_context.append(context);
	m_contextMutex.unlock();

	addRecentFile(context->fileName);

	ui->tabBar->addTab(context->shortName);

	/* Select the newly loaded scene */
	int index = ui->tabBar->count()-1;
	if (ui->tabBar->currentIndex() != index)
		ui->tabBar->setCurrentIndex(index);
	updateUI();
	adjustSize();
}

void MainWindow::addRecentFile(QString fileName) {
	/* Update recent files list */
	QSettings settings("mitsuba-renderer.org", "qtgui");
	QStringList files = settings.value("recentFileList").toStringList();
	files.removeAll(fileName);
	files.prepend(fileName);

	while (files.size() > MAX_RECENT_FILES)
		files.removeLast();

	settings.setValue("recentFileList", files);
	foreach (QWidget *widget, QApplication::topLevelWidgets()) {
		MainWindow *mainWin = qobject_cast<MainWindow *>(widget);
		if (mainWin)
			mainWin->updateRecentFileActions();
	}
}

void MainWindow::updateRecentFileActions() {
	QSettings settings("mitsuba-renderer.org", "qtgui");
	QStringList files = settings.value("recentFileList").toStringList();

	int numRecentFiles = qMin(files.size(), MAX_RECENT_FILES);

	for (int i = 0; i < numRecentFiles; ++i) {
		m_actRecent[i]->setText(QFileInfo(files[i]).fileName());
		m_actRecent[i]->setData(files[i]);
		m_actRecent[i]->setVisible(true);
	}

	for (int j = numRecentFiles; j < MAX_RECENT_FILES; ++j)
		m_actRecent[j]->setVisible(false);
}

void MainWindow::updateUI() {
	int index = ui->tabBar->currentIndex();
	bool hasTab = (index != -1);

	SceneContext *context = hasTab ? m_context[index] : NULL;
	bool isRendering = hasTab ? context->renderJob != NULL : false;
	bool isShowingRendering = hasTab ? context->mode == ERender : false;
	bool hasScene = hasTab && context->scene != NULL;
	bool isInactiveScene = (hasTab && hasScene) ? context->renderJob == NULL : false;
	bool fallback = ui->glView->isUsingSoftwareFallback();

	ui->actionStop->setEnabled(isShowingRendering);
	if (isShowingRendering && !isRendering)
		ui->actionStop->setToolTip(tr("Return to the realtime preview"));
	else
		ui->actionStop->setToolTip(tr("Stop rendering"));

	ui->actionRender->setEnabled(isInactiveScene);
	ui->actionRefresh->setEnabled(isInactiveScene);
	ui->actionRenderSettings->setEnabled(isInactiveScene);
	ui->actionSave->setEnabled(hasScene);
	ui->actionSaveAs->setEnabled(hasScene);
	ui->actionExportImage->setEnabled(hasTab);
	ui->actionClose->setEnabled(hasTab);
	ui->actionDuplicateTab->setEnabled(hasTab);
	ui->actionAdjustSize->setEnabled(hasTab);
#if !defined(__OSX__)
	ui->actionPreviewSettings->setEnabled(!fallback && hasTab);
#else
	bool isVisible = m_previewSettings != NULL && m_previewSettings->isVisible();
	ui->actionPreviewSettings->setEnabled(hasTab && !isVisible && !fallback);
#endif

	if (isRendering) {
		if (!m_progress->isVisible()) {
			QGridLayout *centralLayout = static_cast<QGridLayout *>(centralWidget()->layout());
			centralLayout->addWidget(m_progressWidget, 3, 0, 1, 3);
			m_progressWidget->show();
		}
		m_progress->setValue(context->progress);
		QString etaString;

		if (context->eta != "")
			etaString = QString(" (ETA: ") + context->eta  + QString(")");

		if (context->progressName != "")
			m_progressLabel->setText(QString("<b>%1</b> %2%%3")
				.arg(context->progressName)
				.arg(context->progress, 0, 'f', 1)
				.arg(etaString));
		else
			m_progressLabel->setText("");
	} else {
		if (m_progress->isVisible()) {
			QLayout *centralLayout = centralWidget()->layout();
			m_progressWidget->hide();
			centralLayout->removeWidget(m_progressWidget);
		}
	}
	centralWidget()->updateGeometry();
	layout()->activate();
}

void MainWindow::resizeEvent(QResizeEvent *event) {
	QMainWindow::resizeEvent(event);
	ui->glView->updateScrollBars();
}
	
void MainWindow::on_tabBar_customContextMenuRequested(const QPoint &pt) {
	if (pt.isNull())
		return;
	int tabIndex = ui->tabBar->tabAt(pt);
	if (tabIndex == -1)
		return;
	m_contextIndex = tabIndex;
	QMenu menu(this);
	menu.addAction(ui->actionDuplicateTab);
	if (tabIndex == ui->tabBar->currentIndex())
		menu.addAction(ui->actionAdjustSize);
	menu.addAction(ui->actionClose);
	menu.exec(ui->tabBar->mapToGlobal(pt));
	m_contextIndex = -1;
}

void MainWindow::on_tabBar_tabMoved(int from, int to) {
	m_contextMutex.lock();
	m_context.move(from, to);
	m_contextMutex.unlock();
}

void MainWindow::on_tabBar_currentChanged(int index) {
	if (m_lastTab != NULL) 
		m_lastTab->windowSize = size();

	ui->glView->ignoreResizeEvents(true);
	if (ui->tabBar->currentIndex() != -1)
		ui->tabBar->show();
	else
		ui->tabBar->hide();
	ui->glView->ignoreResizeEvents(false);

	if (index != -1) {
		if (m_context[index] != m_lastTab) {
			ui->glView->setScene(m_context[index]);
#if defined(__OSX__)
			if (m_previewSettings)
				m_previewSettings->setContext(m_context[index]);
#endif
		}
	} else {
		ui->glView->setScene(NULL);
#if defined(__OSX__)
		if (m_previewSettings && m_previewSettings->isVisible())
			m_previewSettings->hide();
#endif
	}
	m_statusMessage = "";
	updateStatus();
	updateUI();

	if (index != -1) {
		const QSize &windowSize = m_context[index]->windowSize;
		if (windowSize.isValid()) {
#if defined(__LINUX__)
			int error = (sizeHint()-windowSize).height();
			if (error > 0 && error <= 5)
				resize(windowSize + QSize(0, error));
			else
				resize(windowSize);
#else
			resize(windowSize);
#endif
		} else {
			adjustSize();
		}

		m_lastTab = m_context[index];
	} else {
		adjustSize();
		m_lastTab = NULL;
	}
	ui->glView->setFocus();
}

void MainWindow::on_actionClose_triggered() {
	int index = ui->tabBar->currentIndex();
	if (m_contextIndex != -1)
		index = m_contextIndex;
	on_tabBar_tabCloseRequested(index);
}

void MainWindow::on_actionAdjustSize_triggered() {
	adjustSize();
}

bool MainWindow::on_tabBar_tabCloseRequested(int index) {
	SceneContext *context = m_context[index];
	if (context->renderJob != NULL) {
		QMessageBox box(QMessageBox::Question, tr("Really close?"), 
			tr("Rendering of scene \"%1\" is unfinished - all progress will be lost. Are you sure you want to continue?").arg(context->shortName),
			QMessageBox::Yes | QMessageBox::No, this, Qt::Dialog | Qt::MSWindowsFixedSizeDialogHint | Qt::Sheet);
		box.setWindowModality(Qt::WindowModal);
		if (box.exec() == QMessageBox::No)
			return false;
	}
	m_contextMutex.lock();
	m_context.removeAt(index);
	m_contextMutex.unlock();
	ui->tabBar->removeTab(index);
	if (context->renderJob != NULL) {
		context->renderJob->cancel();
		context->renderJob->join();
	}
	ui->glView->makeCurrent();
	delete context;
	return true;
}

void MainWindow::closeEvent(QCloseEvent *event) {
	QSettings settings("mitsuba-renderer.org", "qtgui");
	settings.setValue("pos", pos());
	for (int i=ui->tabBar->count()-1; i>=0; --i) {
		if (!on_tabBar_tabCloseRequested(i)) {
			event->ignore();
			return;
		}
	}
	QMainWindow::closeEvent(event);
	m_logWidget->hide();
	Logger *logger = Thread::getThread()->getLogger();
	logger->removeAppender(m_consoleAppender);
	event->accept();
}

bool MainWindow::isActive() {
	if (isActiveWindow() || m_activeWindowHack)
		return true;
	else if (m_currentChild != NULL && m_currentChild->isActiveWindow())
		return true;
#if defined(__OSX__)
	if (m_previewSettings != NULL && m_previewSettings->isActiveWindow()) 
		return true;
#endif
	return false;
}

void MainWindow::drawHLine(SceneContext *ctx, int x1, int y, int x2, const float *color) {
	float *framebuffer = ctx->framebuffer->getFloatData();
	int fbOffset = (x1 + y*ctx->framebuffer->getWidth())*4;
	for (int x=x1; x<=x2; x++) {
		framebuffer[fbOffset] = color[0];
		framebuffer[fbOffset+1] = color[1];
		framebuffer[fbOffset+2] = color[2];
		fbOffset+=4;
	}
}

void MainWindow::drawVLine(SceneContext *ctx, int x, int y1, int y2, const float *color) {
	float *framebuffer = ctx->framebuffer->getFloatData();
	int width = ctx->framebuffer->getWidth(), fbOffset = (x + y1*width)*4;
	for (int y=y1; y<=y2; y++) {
		framebuffer[fbOffset] = color[0];
		framebuffer[fbOffset+1] = color[1];
		framebuffer[fbOffset+2] = color[2];
		fbOffset += width * 4;
	}
}

void MainWindow::on_actionRenderSettings_triggered() {
	int currentIndex = ui->tabBar->currentIndex();
	SceneContext *context = m_context[currentIndex];
	RenderSettingsDialog *dialog = new RenderSettingsDialog(this);
	dialog->setAttribute(Qt::WA_DeleteOnClose);
	dialog->setWindowModality(Qt::WindowModal);
	dialog->load(context);
	connect(dialog, SIGNAL(finished(int)), this, SLOT(onRenderSettingsClose(int)));
	m_currentChild = dialog;
	// prevent a tab drawing artifact on Qt/OSX
	m_activeWindowHack = true;
	dialog->show();
	qApp->processEvents();
	m_activeWindowHack = false;
}

void MainWindow::onRenderSettingsClose(int reason) {
	RenderSettingsDialog *dialog = static_cast<RenderSettingsDialog *>(sender());
	int currentIndex = ui->tabBar->currentIndex();
	SceneContext *context = m_context[currentIndex];
	m_currentChild = NULL;
	if (reason == QDialog::Accepted) {
		Scene *scene = context->scene;
		scene->incRef();
		dialog->apply(context);
		ui->glView->setPathLength(context->detectPathLength());
		if (dialog->resolutionHasChanged()) {
			ui->glView->refreshScene();
			updateUI();
			adjustSize();
		}
		scene->decRef();
#if defined(__OSX__)
		if (m_previewSettings)
			m_previewSettings->setContext(context);
#endif
	}
}

void MainWindow::on_actionPreviewSettings_triggered() {
#if !defined(__OSX__)
	SceneContext *context = m_context[ui->tabBar->currentIndex()];
	PreviewSettingsDialog d(this, context, ui->glView->getRendererCapabilities());
	connect(&d, SIGNAL(pathLengthChanged(int)), ui->glView, SLOT(setPathLength(int)));
	connect(&d, SIGNAL(clampingChanged(Float)), ui->glView, SLOT(setClamping(Float)));
	connect(&d, SIGNAL(shadowMapResolutionChanged(int)), ui->glView, SLOT(setShadowMapResolution(int)));
	connect(&d, SIGNAL(gammaChanged(bool, Float)), ui->glView, SLOT(setGamma(bool, Float)));
	connect(&d, SIGNAL(exposureChanged(Float)), ui->glView, SLOT(setExposure(Float)));
	connect(&d, SIGNAL(reinhardKeyChanged(Float)), ui->glView, SLOT(setReinhardKey(Float)));
	connect(&d, SIGNAL(reinhardBurnChanged(Float)), ui->glView, SLOT(setReinhardBurn(Float)));
	connect(&d, SIGNAL(previewMethodChanged(EPreviewMethod)), ui->glView, SLOT(setPreviewMethod(EPreviewMethod)));
	connect(&d, SIGNAL(toneMappingMethodChanged(EToneMappingMethod)), ui->glView, SLOT(setToneMappingMethod(EToneMappingMethod)));
	connect(&d, SIGNAL(diffuseReceiversChanged(bool)), ui->glView, SLOT(setDiffuseReceivers(bool)));
	connect(&d, SIGNAL(diffuseSourcesChanged(bool)), ui->glView, SLOT(setDiffuseSources(bool)));
	d.setMaximumSize(d.minimumSize());
	d.exec();
	QSettings settings("mitsuba-renderer.org", "qtgui");
	settings.setValue("preview_sRGB", context->srgb);
	settings.setValue("preview_gamma", context->gamma);
	settings.setValue("preview_exposure", context->exposure);
	settings.setValue("preview_shadowMapResolution", context->shadowMapResolution);
	settings.setValue("preview_reinhardKey", context->reinhardKey);
	settings.setValue("preview_reinhardBurn", context->reinhardBurn);
	settings.setValue("preview_clamping", context->clamping);
	settings.setValue("preview_previewMethod", context->previewMethod);
	settings.setValue("preview_toneMappingMethod", context->toneMappingMethod);
	settings.setValue("preview_diffuseReceivers", context->diffuseReceivers);
	settings.setValue("preview_diffuseSources", context->diffuseSources);
#else
	if (!m_previewSettings) {
		m_previewSettings = new PreviewSettingsDlg(this);
		connect(m_previewSettings, SIGNAL(pathLengthChanged(int)), ui->glView, SLOT(setPathLength(int)));
		connect(m_previewSettings, SIGNAL(clampingChanged(Float)), ui->glView, SLOT(setClamping(Float)));
		connect(m_previewSettings, SIGNAL(shadowMapResolutionChanged(int)), ui->glView, SLOT(setShadowMapResolution(int)));
		connect(m_previewSettings, SIGNAL(gammaChanged(bool, Float)), ui->glView, SLOT(setGamma(bool, Float)));
		connect(m_previewSettings, SIGNAL(exposureChanged(Float)), ui->glView, SLOT(setExposure(Float)));
		connect(m_previewSettings, SIGNAL(reinhardKeyChanged(Float)), ui->glView, SLOT(setReinhardKey(Float)));
		connect(m_previewSettings, SIGNAL(reinhardBurnChanged(Float)), ui->glView, SLOT(setReinhardBurn(Float)));
		connect(m_previewSettings, SIGNAL(previewMethodChanged(EPreviewMethod)), ui->glView, SLOT(setPreviewMethod(EPreviewMethod)));
		connect(m_previewSettings, SIGNAL(toneMappingMethodChanged(EToneMappingMethod)), ui->glView, SLOT(setToneMappingMethod(EToneMappingMethod)));
		connect(m_previewSettings, SIGNAL(close()), this, SLOT(onPreviewSettingsClose()));
		connect(m_previewSettings, SIGNAL(diffuseReceiversChanged(bool)), ui->glView, SLOT(setDiffuseReceivers(bool)));
		connect(m_previewSettings, SIGNAL(diffuseSourcesChanged(bool)), ui->glView, SLOT(setDiffuseSources(bool)));
	}
	SceneContext *ctx = NULL;
	if (ui->tabBar->currentIndex() != -1)
		ctx = m_context[ui->tabBar->currentIndex()];
	m_previewSettings->setContext(ctx);
	m_previewSettings->show();
	ui->actionPreviewSettings->setEnabled(false);
#endif
}

void MainWindow::onPreviewSettingsClose() {
	int index = ui->tabBar->currentIndex();
	bool hasTab = index != -1;
	ui->actionPreviewSettings->setEnabled(hasTab);
	if (hasTab) {
		SceneContext *context = m_context[index];
		QSettings settings("mitsuba-renderer.org", "qtgui");
		settings.setValue("preview_sRGB", context->srgb);
		settings.setValue("preview_gamma", context->gamma);
		settings.setValue("preview_exposure", context->exposure);
		settings.setValue("preview_shadowMapResolution", context->shadowMapResolution);
		settings.setValue("preview_reinhardKey", context->reinhardKey);
		settings.setValue("preview_reinhardBurn", context->reinhardBurn);
		settings.setValue("preview_clamping", context->clamping);
		settings.setValue("preview_previewMethod", context->previewMethod);
		settings.setValue("preview_toneMappingMethod", context->toneMappingMethod);
		settings.setValue("preview_diffuseSources", context->diffuseSources);
		settings.setValue("preview_diffuseReceivers", context->diffuseReceivers);
	}
}

void MainWindow::on_actionSettings_triggered() {
	Logger *logger = Thread::getThread()->getLogger();
	ref<Scheduler> sched = Scheduler::getInstance();
	std::vector<Worker *> localWorkers;
	
	if (m_renderQueue->getJobCount() != 0) {
		QMessageBox::warning(this, tr("Rendering in progress"),
			tr("The program settings cannot be changed while a rendering is in progress."),
			QMessageBox::Ok);
		return;
	}

	int workerCount = sched->getWorkerCount();
	for (int i=0; i<workerCount; ++i) {
		Worker *worker = sched->getWorker(i);
		if (worker->getClass()->derivesFrom(LocalWorker::m_theClass))
			localWorkers.push_back(worker);
	}

	ProgramSettingsDialog d(this);
	d.setWindowModality(Qt::ApplicationModal);
	d.setLogLevel(logger->getLogLevel());
	d.setInvertMouse(ui->glView->getInvertMouse());
	d.setNavigationMode(ui->glView->getNavigationMode());
	d.setMouseSensitivity(ui->glView->getMouseSensitivity());
	d.setBlockSize(m_blockSize);
	d.setSearchPaths(m_searchPaths);
	d.setLocalWorkerCount(localWorkers.size());
	d.setConnections(m_connections);
	d.setCheckForUpdates(m_checkForUpdates);
	d.setListenPort(m_listenPort);
	d.setNodeName(m_nodeName);

	if (d.exec()) {
		QList<ServerConnection> &connections = d.getConnections();
		QList<QVariant> connectionData;
		for (int i=0; i<connections.size(); ++i)
			connectionData.append(connections[i].toByteArray());

		QSettings settings("mitsuba-renderer.org", "qtgui");
		settings.setValue("verbosity", d.getLogLevel());
		settings.setValue("invertMouse", d.getInvertMouse());
		settings.setValue("blockSize", d.getBlockSize());
		settings.setValue("searchPaths", d.getSearchPaths());
		settings.setValue("localWorkers", d.getLocalWorkerCount());
		settings.setValue("connections", connectionData);
		settings.setValue("checkForUpdates", d.getCheckForUpdates());
		settings.setValue("mouseSensitivity", d.getMouseSensitivity());
		settings.setValue("listenPort", d.getListenPort());
		settings.setValue("nodeName", d.getNodeName());
		settings.setValue("navigationMode", (int) d.getNavigationMode());

		logger->setLogLevel(d.getLogLevel());
		ui->glView->setInvertMouse(d.getInvertMouse());
		ui->glView->setMouseSensitivity(d.getMouseSensitivity());
		ui->glView->setNavigationMode(d.getNavigationMode());
		m_blockSize = d.getBlockSize();
		m_searchPaths = d.getSearchPaths();
		m_checkForUpdates = d.getCheckForUpdates();
		m_listenPort = d.getListenPort();
		m_nodeName = d.getNodeName();

		bool localWorkersChanged = (int) localWorkers.size() != d.getLocalWorkerCount();

		if (localWorkersChanged || m_connections != d.getConnections()) {
			ref<Scheduler> sched = Scheduler::getInstance();
			sched->pause();
			while (d.getLocalWorkerCount() > (int) localWorkers.size()) {
				LocalWorker *worker = new LocalWorker(formatString("wrk%i", localWorkerCtr++));
				sched->registerWorker(worker);
				localWorkers.push_back(worker);
			}
			while (d.getLocalWorkerCount() < (int) localWorkers.size()) {
				Worker *worker = localWorkers.back();
				sched->unregisterWorker(worker);
				localWorkers.pop_back();
			}
			QList<ServerConnection> removeList, 
				&newConnections = d.getConnections();
			for (int i=0; i<m_connections.size(); ++i) {
				ServerConnection &c = m_connections[i];
				if (!newConnections.contains(c)) 
					removeList.append(c);
			}
			for (int i=0; i<newConnections.size(); ++i) {
				ServerConnection &c = newConnections[i];
				if (!m_connections.contains(c)) {
					sched->registerWorker(c.worker);
					c.isRegistered = true;
					c.worker->decRef();
					m_connections.append(c);
				}
			}
			for (int i=0; i<removeList.size(); ++i) {
				ServerConnection &c = removeList[i];
				sched->unregisterWorker(c.worker);
				m_connections.removeAll(c);
			}
			sched->start();
		}
	}
}


void MainWindow::on_actionStop_triggered() {
	SceneContext *context = m_context[ui->tabBar->currentIndex()];
	if (context->renderJob) {
		context->cancelled = true;
		context->renderJob->cancel();
	} else if (context->mode == ERender) {
		context->mode = EPreview;
	}
	updateUI();
}

void MainWindow::on_actionShowLog_triggered() {
	m_logWidget->show();
}

void MainWindow::on_actionRender_triggered() {
	int index = ui->tabBar->currentIndex();
	if (index == -1)
		return;
	SceneContext *context = m_context[index];
	if (context->renderJob != NULL)
		return;

	Scene *scene = context->scene;
	scene->setBlockSize(m_blockSize);
	context->renderJob = new RenderJob("rend", scene, m_renderQueue, NULL, 
		context->sceneResID, -1, -1, false);
	context->cancelMode = ERender;
	if (context->mode != ERender)
		ui->glView->downloadFramebuffer();
	context->cancelled = false;
	context->progress = 0;
	context->eta = "";
	context->progressName = "";
	context->mode = ERender;
	m_statusMessage = "";
	updateUI();
	context->sizeIncrease = QSize(0, m_progressWidget->sizeHint().height());

#if defined(__LINUX__)
	/* Workaround: on Linux, a few pixels get lost somehow, which
	   otherwise causes an unnecessary scrollbar to appear */
	int error = (sizeHint()-size()).height() - context->sizeIncrease.height();
	if (error > 0 && error <= 5)
		resize(size() + context->sizeIncrease + QSize(0, error));
	else
		resize(size() + context->sizeIncrease);
#else
	resize(size() + context->sizeIncrease);
#endif

	updateStatus();
	context->renderJob->start();
}

void MainWindow::on_actionRefresh_triggered() {
	int index = ui->tabBar->currentIndex();
	if (index < 0)
		return;
	SceneContext *context = m_context[index];
	on_tabBar_currentChanged(-1);
	SceneContext *newContext = loadScene(context->fileName);
	if (newContext == NULL) 
		return;
	delete context;
	m_contextMutex.lock();
	m_context[index] = newContext;
	m_contextMutex.unlock();
	qApp->processEvents();
	on_tabBar_currentChanged(index);
}

inline float toSRGB(float value) {
	if (value < 0.0031308f)
		return 12.92f * value;
	return 1.055f * std::pow(value, 0.41666f) - 0.055f;
}

void MainWindow::on_actionExportImage_triggered() {
	SceneContext *ctx = m_context[ui->tabBar->currentIndex()];
	QFileDialog dialog(this, tr("Export image .."),
		"", tr("Linear EXR Image (*.exr);; Tonemapped 8-bit PNG Image (*.png)"));

	QSettings settings("mitsuba-renderer.org", "qtgui");
	dialog.setViewMode(QFileDialog::Detail);
	dialog.setAcceptMode(QFileDialog::AcceptSave);

#if defined(__OSX__)
	dialog.setOption(QFileDialog::DontUseNativeDialog, true);
#endif

	dialog.restoreState(settings.value("fileDialogState").toByteArray());

    if (dialog.exec() == QDialog::Accepted) {
        QString fileName = dialog.selectedFiles().value(0);
		Bitmap::EFileFormat format;
		settings.setValue("fileDialogState", dialog.saveState());

		if (fileName.endsWith(".exr")) {
			format = Bitmap::EEXR;
		} else if (fileName.endsWith(".png")) {
			format = Bitmap::EPNG;
		} else {
			SLog(EError, "Unknown file type -- the filename must end in either .exr or .png");
			return;
		}

		ref<FileStream> fs = new FileStream(qPrintable(fileName), FileStream::ETruncReadWrite);

		if (ctx->mode == EPreview)
			ui->glView->downloadFramebuffer();

		if (format == Bitmap::EEXR) {
			ctx->framebuffer->save(format, fs);
		} else {
			int width = ctx->framebuffer->getWidth();
			int height = ctx->framebuffer->getHeight();
			ref<Bitmap> temp = new Bitmap(width, height, 32);
			float *source = ctx->framebuffer->getFloatData();
			uint8_t *target = temp->getData();
			float invGamma = 1.0f/(float) ctx->gamma;
			float exposure = std::pow(2.0f, (float) ctx->exposure);
			Float reinhardKey = 0;
			Float invWpSqr = std::pow((Float) 2, (Float) ctx->reinhardBurn);
			if (ctx->toneMappingMethod == EReinhard) {
				Float avgLogLuminance = 0;
				for (int y=0; y<height; ++y) {
					for (int x=0; x<width; ++x) {
						Spectrum spec;
						spec.fromLinearRGB(source[(y*width+x)*4+0], source[(y*width+x)*4+1], source[(y*width+x)*4+2]);
						avgLogLuminance += std::log(0.001f+spec.getLuminance());
					}
				}
				avgLogLuminance = std::exp(avgLogLuminance/(width*height));
				reinhardKey = ctx->reinhardKey / avgLogLuminance;
			}

			for (int y=0; y<height; ++y) {
				for (int x=0; x<width; ++x) {
					float r, g, b, a = source[3];

					if (ctx->toneMappingMethod == EGamma) {
						r = source[0]*exposure;
						g = source[1]*exposure;
						b = source[2]*exposure;
					} else {
						Spectrum spec;
						Float X, Y, Z;
						spec.fromLinearRGB(source[0], source[1], source[2]);
						spec.toXYZ(X, Y, Z);
						Float normalization = 1/(X + Y + Z);
						Float x = X*normalization, y = Y*normalization;
						Float Lp = Y * reinhardKey;
						Y = Lp * (1.0f + Lp*invWpSqr) / (1.0f + Lp);
						X = x * (Y/y); 
						Z = (Y/y) * (1.0 - x - y);
						spec.fromXYZ(X, Y, Z);
						Float rF, gF, bF;
						spec.toLinearRGB(rF, gF, bF);
						r = rF; g = gF; b = bF; 
					}

					if (ctx->srgb) {
						r = toSRGB(r);
						g = toSRGB(g);
						b = toSRGB(b);
					} else {
						r = std::pow(r, invGamma);
						g = std::pow(g, invGamma);
						b = std::pow(b, invGamma);
					}

					*target++ = (uint8_t) std::min(255, std::max(0, (int) (r*255)));
					*target++ = (uint8_t) std::min(255, std::max(0, (int) (g*255)));
					*target++ = (uint8_t) std::min(255, std::max(0, (int) (b*255)));
					*target++ = (uint8_t) std::min(255, std::max(0, (int) (a*255)));
					source += 4;
				}
			}
			temp->setGamma(ctx->srgb ? -1 : ctx->gamma);
			temp->save(format, fs);
		}
	}
}

void MainWindow::on_actionSave_triggered() {
	SceneContext *context = m_context[ui->tabBar->currentIndex()];
	saveScene(this, context, context->fileName);
}

void MainWindow::on_actionSaveAs_triggered() {
	QFileDialog *dialog = new QFileDialog(this, tr("Save as .."),
		"", tr("Mitsuba scenes (*.xml)"));

	m_currentChild = dialog;
	QSettings settings("mitsuba-renderer.org", "qtgui");
	dialog->setAttribute(Qt::WA_DeleteOnClose);
	dialog->setViewMode(QFileDialog::Detail);
	dialog->setAcceptMode(QFileDialog::AcceptSave);
	dialog->restoreState(settings.value("fileDialogState").toByteArray());
	dialog->setWindowModality(Qt::WindowModal);
	connect(dialog, SIGNAL(finished(int)), this, SLOT(onSaveAsDialogClose(int)));
	m_currentChild = dialog;
	// prevent a tab drawing artifact on Qt/OSX
	m_activeWindowHack = true;
	dialog->show();
	qApp->processEvents();
	m_activeWindowHack = false;
}

void MainWindow::onSaveAsDialogClose(int reason) {
	int currentIndex = ui->tabBar->currentIndex();
	SceneContext *context = m_context[currentIndex];

	QSettings settings("mitsuba-renderer.org", "qtgui");
	QFileDialog *dialog = static_cast<QFileDialog *>(sender());
	m_currentChild = NULL;
	if (reason == QDialog::Accepted) {
		m_currentChild = NULL;
        QString fileName = dialog->selectedFiles().value(0);
		settings.setValue("fileDialogState", dialog->saveState());
		saveScene(this, context, fileName);
		context->fileName = fileName;
		context->shortName = QFileInfo(fileName).fileName();
		ui->tabBar->setTabText(currentIndex, context->shortName);
		addRecentFile(fileName);
	}
}

void MainWindow::on_actionNavigationControls_triggered() {
	NavigationDialog nav(this);
	nav.exec();
}

void MainWindow::on_actionReferenceManual_triggered() {
	QDesktopServices::openUrl(QUrl("http://www.mitsuba-renderer.org/docs.html"));
}
	
void MainWindow::on_actionAbout_triggered() {
	AboutDialog about(this);
	about.exec();
}

void MainWindow::onJobFinished(const RenderJob *job, bool cancelled) {
	SceneContext *context = getContext(job, false);
	m_renderQueue->join();
	if (context == NULL)
		return;
	if (cancelled) {
		if (!context->cancelled) {
			QMessageBox::critical(this, tr("Error while rendering"),
				tr("The rendering job did not complete successfully. Please check the log."), 
				QMessageBox::Ok);
		} else {
			context->mode = context->cancelMode;
			if (ui->tabBar->currentIndex() != -1 &&
				m_context[ui->tabBar->currentIndex()] == context)
				ui->glView->resumePreview();
		}
	}
	onRefresh(job);
	context->renderJob = NULL;
	updateUI();
	if (ui->tabBar->currentIndex() != -1 &&
		m_context[ui->tabBar->currentIndex()] == context)
		resize(size() - context->sizeIncrease);
}
	
void MainWindow::onStatusMessage(const QString &status) {
	m_statusMessage = status;
	updateStatus();
}

void MainWindow::updateStatus() {
	if (m_statusMessage == "")
		setWindowTitle(tr("Mitsuba renderer"));
	else
		setWindowTitle(tr("Mitsuba renderer [%1]").arg(m_statusMessage));
}
    
void MainWindow::on_actionStartServer_triggered() {
	m_serverWidget = new ServerWidget(NULL, m_nodeName, m_listenPort);
	ui->actionStartServer->setEnabled(false);
	connect(m_serverWidget, SIGNAL(closed()), this, SLOT(onServerClosed()));
	m_serverWidget->show();
}

void MainWindow::onServerClosed() {
	delete m_serverWidget;
	m_serverWidget = NULL;
	ui->actionStartServer->setEnabled(true);
}

void MainWindow::onWorkBegin(const RenderJob *job, const RectangularWorkUnit *wu, int worker) {
	SceneContext *context = getContext(job, false);
	if (context == NULL)
		return;

	Film *film = context->scene->getFilm();
	int ox = wu->getOffset().x - film->getCropOffset().x, 
		oy = wu->getOffset().y - film->getCropOffset().y,
		ex = ox + wu->getSize().x, ey = oy + wu->getSize().y;
	const float *color = NULL;

	/* Use desaturated colors to highlight the host
	   responsible for rendering the current image block */
	const float white[]     = { 1.0f, 1.0f, 1.0f };
	const float red[]       = { 1.0f, 0.3f, 0.3f };
	const float green[]     = { 0.3f, 1.0f, 0.3f };
	const float blue[]      = { 0.3f, 0.3f, 1.0f };
	const float gray[]      = { 0.5f, 0.5f, 0.5f };
	const float yellow[]    = { 1.0f, 1.0f, 0.0f };
	const float magenta[]   = { 1.0f, 0.3f, 1.0f };
	const float turquoise[] = { 0.3f, 1.0f, 1.0f };

	switch (worker % 8) {
		case 1: color = green; break;
		case 2: color = yellow; break;
		case 3: color = blue; break;
		case 4: color = gray; break;
		case 5: color = red; break;
		case 6: color = magenta; break;
		case 7: color = turquoise; break;
		case 0:
		default:
			color = white;
			break;
	}

	if (wu->getSize().x < 3 || wu->getSize().y < 3)
		return;

	drawHLine(context, ox, oy, ox + 3, color);
	drawHLine(context, ex - 4, oy, ex - 1, color);
	drawHLine(context, ox, ey - 1, ox + 3, color);
	drawHLine(context, ex - 4, ey - 1, ex - 1, color);
	drawVLine(context, ox, oy, oy + 3, color);
	drawVLine(context, ex - 1, oy, oy + 3, color);
	drawVLine(context, ex - 1, ey - 4, ey - 1, color);
	drawVLine(context, ox, ey - 4, ey - 1, color);

	/* This is not executed in the event loop -- take some precautions */
	m_contextMutex.lock();
	bool isCurrentView = ui->tabBar->currentIndex() < m_context.size() &&
		m_context[ui->tabBar->currentIndex()] == context;
	m_contextMutex.unlock();
	if (isCurrentView)
		emit updateView();
}
	
void MainWindow::on_glView_loadFileRequest(const QString &string) {
	loadFile(string);
}

void MainWindow::onWorkEnd(const RenderJob *job, const ImageBlock *block) {
	int ox = block->getOffset().x, oy = block->getOffset().y,
		ex = ox + block->getSize().x, ey = oy + block->getSize().y;
	SceneContext *context = getContext(job, false);
	if (context == NULL)
		return;
	Film *film = context->scene->getFilm();
	Point2i co = film->getCropOffset();
	Bitmap *bitmap = context->framebuffer;
	float *framebuffer = bitmap->getFloatData();
	Float r, g, b;

	for (int y = oy; y < ey; ++y) {
		int fbOffset = (ox - co.x + (y - co.y)*bitmap->getWidth())*4;
		for (int x = ox; x < ex; ++x) {
			film->getValue(x, y).toLinearRGB(r, g, b);
			framebuffer[fbOffset] = (float) r;
			framebuffer[fbOffset+1] = (float) g;
			framebuffer[fbOffset+2] = (float) b;
			framebuffer[fbOffset+3] = 1;
			fbOffset += 4;
		}
	}

	/* This is executed by worker threads -- take some precautions */
	m_contextMutex.lock();
	bool isCurrentView = ui->tabBar->currentIndex() < m_context.size() &&
		m_context[ui->tabBar->currentIndex()] == context;
	m_contextMutex.unlock();
	if (isCurrentView)
		emit updateView();
}

void MainWindow::onRefresh(const RenderJob *job) {
	SceneContext *context = getContext(job, false);
	if (context == NULL)
		return;
	Film *film = context->scene->getFilm();
	Point2i co = film->getCropOffset();
	Bitmap *bitmap = context->framebuffer;
	float *framebuffer = bitmap->getFloatData();
	Float r, g, b;

	for (int y = 0; y < bitmap->getHeight(); ++y) {
		int fbOffset = y*bitmap->getWidth()*4;
		for (int x = 0; x < bitmap->getWidth(); ++x) {
			film->getValue(x + co.x, y + co.y).toLinearRGB(r, g, b);
			framebuffer[fbOffset] = (float) r;
			framebuffer[fbOffset+1] = (float) g;
			framebuffer[fbOffset+2] = (float) b;
			fbOffset += 4;
		}
	}
	
	if (m_context[ui->tabBar->currentIndex()] == context)
		emit updateView();
}

bool ServerConnection::createWorker(QWidget *parent) {
	ref<Stream> stream;
	try {
		if (type == EDirectConnection) {
			stream = new SocketStream(hostName.toStdString(), port);
		} else {
			std::vector<std::string> cmdLine;
			cmdLine.push_back(formatString("bash -c 'cd %s; . setpath.sh; mtssrv -ls'", instDir.toLatin1().constData()));
			stream = new SSHStream(userName.toStdString(), 
				hostName.toStdString(), cmdLine, port);
		}
		worker = new RemoteWorker(formatString("net%i", remoteWorkerCtr++), stream);
		return true;
	} catch (const std::exception &e) {
		QString extra;
		if (type == ESSHConnection && !QString(e.what()).contains("configuration mismatch"))
			extra = parent->tr(" - Please make sure that you can log in manually using the "
				"command line and that <a href=\"http://www.debian-administration.org/articles/152\">"
				"passwordless authentication</a> is active.");
		QMessageBox::critical(parent, parent->tr("Unable to connect to %1").arg(hostName),
			QString("Unable to create a connection to \"%1\": %2%3")
			.arg(hostName).arg(e.what()).arg(extra), QMessageBox::Ok);
		return false;
	}
}

QSize MainWindow::sizeHint() const {
	QSize hint = QMainWindow::sizeHint();
	/* Don't include scroll bars in the size hint */
	if (ui->hScrollBar->isVisible())
		hint -= QSize(0, ui->hScrollBar->sizeHint().height());
	if (ui->vScrollBar->isVisible())
		hint -= QSize(ui->vScrollBar->sizeHint().width(), 0);
	return hint;
}

QString ServerConnection::toString() const {
	return QString("%1 (%2, %3 cores)")
		.arg(worker->getNodeName().c_str())
		.arg(type == EDirectConnection ? "direct" : "ssh")
		.arg(worker->getCoreCount());
}

SceneContext::SceneContext(SceneContext *ctx) {
	if (ctx->scene) {
		scene = new Scene(ctx->scene);
		ref<PluginManager> pluginMgr = PluginManager::getInstance();
		ref<PinholeCamera> oldCamera = static_cast<PinholeCamera *>(ctx->scene->getCamera());
		ref<PinholeCamera> camera = static_cast<PinholeCamera *> 
			(pluginMgr->createObject(Camera::m_theClass, oldCamera->getProperties()));
		ref<Sampler> sampler = static_cast<Sampler *> 
			(pluginMgr->createObject(Sampler::m_theClass, ctx->scene->getSampler()->getProperties()));
		ref<Film> film = static_cast<Film *> 
			(pluginMgr->createObject(Film::m_theClass, oldCamera->getFilm()->getProperties()));
		const Integrator *oldIntegrator = ctx->scene->getIntegrator();
		ref<Integrator> currentIntegrator;

		int depth = 0;
		std::vector<Integrator *> integratorList;
		while (oldIntegrator != NULL) {
			ref<Integrator> integrator = static_cast<Integrator *> (pluginMgr->createObject(
				Integrator::m_theClass, oldIntegrator->getProperties()));
			if (depth++ == 0) 
				scene->setIntegrator(integrator);
			else 
				currentIntegrator->addChild("", integrator);
			currentIntegrator = integrator;
			integratorList.push_back(integrator);
			oldIntegrator = oldIntegrator->getSubIntegrator();
		}

		for (int i=(int) integratorList.size()-1; i>=0; --i)
			integratorList[i]->configure();

		ref<ReconstructionFilter> rfilter = static_cast<ReconstructionFilter *> 
			(pluginMgr->createObject(ReconstructionFilter::m_theClass, oldCamera->getFilm()->
				getReconstructionFilter()->getProperties()));

		rfilter->configure();
		film->addChild("", rfilter);
		film->configure();
		sampler->configure();
		camera->addChild("", sampler);
		camera->addChild("", film);
		camera->setViewTransform(oldCamera->getViewTransform());
		camera->setFov(oldCamera->getFov());
		camera->configure();
		scene->setCamera(camera);
		scene->setSampler(sampler);
		scene->configure();
		sceneResID = ctx->sceneResID;
		Scheduler::getInstance()->retainResource(sceneResID);
	} else {
		sceneResID = -1;
		renderJob = NULL;
	}
	fileName = ctx->fileName;
	shortName = ctx->shortName;
	movementScale = ctx->movementScale;
	up = ctx->up;
	renderJob = NULL;
	cancelled = false;
	progress = 0.0f;
	framebuffer = ctx->framebuffer->clone();
	mode = ctx->renderJob ? EPreview : ctx->mode;
	gamma = ctx->gamma;
	exposure = ctx->exposure;
	clamping = ctx->clamping;
	srgb = ctx->srgb;
	pathLength = ctx->pathLength;
	shadowMapResolution = ctx->shadowMapResolution;
	previewMethod = ctx->previewMethod;
	toneMappingMethod = ctx->toneMappingMethod;
	windowSize = ctx->windowSize;
	sizeIncrease = ctx->sizeIncrease;
	scrollOffset = ctx->scrollOffset;
	reinhardKey = ctx->reinhardKey;
	reinhardBurn = ctx->reinhardBurn;
	diffuseReceivers = ctx->diffuseReceivers;
	diffuseSources = ctx->diffuseSources;
}

SceneContext::~SceneContext() {
	if (scene && sceneResID != -1)
		Scheduler::getInstance()->unregisterResource(sceneResID);
	if (previewBuffer.buffer) {
		previewBuffer.buffer->disassociate();
		previewBuffer.buffer->decRef();
	}
	if (previewBuffer.sync) 
		previewBuffer.sync->decRef();
}

int SceneContext::detectPathLength() const {
	if (!scene)
		return 2;

	const Integrator *integrator = scene->getIntegrator();
	int extraDepth = 0;

	while (integrator->getSubIntegrator() != NULL) {
		if (integrator->getClass()->getName() == "IrradianceCacheIntegrator")
			extraDepth = 1;
		integrator = integrator->getSubIntegrator();
	}
	const Properties &integratorProps = integrator->getProperties();
	int maxDepth = -1;

	if (integratorProps.hasProperty("maxDepth"))
		maxDepth = integratorProps.getInteger("maxDepth");

	if (maxDepth == -1) {
		if (integratorProps.getPluginName() == "direct")
			maxDepth = 2;
		else
			maxDepth = 5;
	}

	return std::max(2, std::min(maxDepth + extraDepth, 6));
}

MTS_IMPLEMENT_CLASS(QRenderListener, false, RenderListener)
