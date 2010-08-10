#include "rendersettingsdlg.h"
#include "ui_rendersettingsdlg.h"
#include <mitsuba/core/plugin.h>

/* ====================== Some helper routines ====================== */

static void setComboBox(QComboBox *box, const std::string &pluginName) {
	for (int i=0; i<box->count(); ++i) {
		const QList<QVariant> &data = box->itemData(i).toList();
		if (data.at(2) == pluginName.c_str()) {
			box->setCurrentIndex(i);
			return;
		}
	}
	SLog(EError, "Unable to find combo box entry named \"%s\"", pluginName.c_str());
}

static std::string getPluginName(QComboBox *box) {
	return box->itemData(box->currentIndex()).toList().at(2).toString().toStdString();
}

/* ====================== RenderSettingsDialog impl ====================== */

RenderSettingsDialog::RenderSettingsDialog(QWidget *parent) :
		QDialog(parent, Qt::Sheet),
	ui(new Ui::RenderSettingsDialog), m_icNode(NULL), m_aiNode(NULL) {
	ui->setupUi(this);

	connect(ui->integratorBox, SIGNAL(highlighted(int)), SLOT(cbHighlighted(int)));
	connect(ui->integratorBox, SIGNAL(activated(int)), SLOT(update()));
	connect(ui->samplerBox, SIGNAL(highlighted(int)), SLOT(cbHighlighted(int)));
	connect(ui->samplerBox, SIGNAL(activated(int)), SLOT(update()));
	connect(ui->rFilterBox, SIGNAL(highlighted(int)), SLOT(cbHighlighted(int)));
	connect(ui->rFilterBox, SIGNAL(activated(int)), SLOT(update()));
	connect(ui->icBox, SIGNAL(pressed()), SLOT(chkBoxPressed()));
	connect(ui->aiBox, SIGNAL(pressed()), SLOT(chkBoxPressed()));
	connect(ui->icBox, SIGNAL(toggled(bool)), SLOT(update()));
	connect(ui->aiBox, SIGNAL(toggled(bool)), SLOT(update()));
	connect(ui->resolutionBox, SIGNAL(activated(int)), SLOT(refresh()));
	connect(ui->resolutionBox, SIGNAL(editTextChanged(const QString &)), SLOT(refresh()));

	QFile file(":/resources/docs.xml");
	if (!file.open(QIODevice::ReadOnly) || !m_document.setContent(&file))
		SLog(EError, "Unable to read the documentation file!");
	file.close();

	/* Populate the integrator, rec. filter & sampler combo box widgets */
	QDomElement docRoot = m_document.documentElement();
	for (QDomElement e = docRoot.firstChildElement("plugin"); !e.isNull();
		 e = e.nextSiblingElement("plugin")) {
		QString docString, name = e.attribute("name");
		if (!e.firstChildElement("descr").isNull()) {
			/* Create a HTML-based documentation string */
			QDomDocument helpDoc;
			QDomElement root = helpDoc.createElement("p");
			helpDoc.appendChild(root);

			for (QDomNode child = e.firstChildElement("descr").firstChild();
			   !child.isNull(); child = child.nextSibling())
			root.appendChild(helpDoc.importNode(child, true));
			docString = helpDoc.toString();
		}

		if (e.hasAttribute("show") && e.attribute("show") == "true") {
			QString type = e.attribute("type"),
					className = e.attribute("className"),
					readableName = e.attribute("readableName"),
					name = e.attribute("name");

			QList<QVariant> list;
			list.append(className);
			list.append(docString);
			list.append(name);

			if (type == "integrator")
				ui->integratorBox->addItem(readableName, list);
			else if (type == "sampler")
				ui->samplerBox->addItem(readableName, list);
			else if (type == "rfilter")
				ui->rFilterBox->addItem(readableName, list);
		}
		if (name == "irrcache")
			ui->icBox->setProperty("help", docString);
		else if (name == "errctrl")
			ui->aiBox->setProperty("help", docString);
	}

	m_model = new XMLTreeModel(docRoot, palette(), this);
	ui->treeView->setModel(m_model);
	ui->treeView->setAlternatingRowColors(true);
	ui->treeView->setUniformRowHeights(true);
	ui->treeView->setColumnWidth(0, 270);
	ui->treeView->setItemDelegate(new PropertyDelegate(this));
	connect(ui->treeView->selectionModel(), SIGNAL(selectionChanged(const QItemSelection &, const QItemSelection)), 
		SLOT(onTreeSelectionChange(const QItemSelection &, const QItemSelection &)));
	connect(m_model, SIGNAL(dataChanged(const QModelIndex&, const QModelIndex&)), this, SLOT(dataChanged()));
	m_integratorNode = m_model->registerClass("MIPathTracer", "Path tracer");
	m_samplerNode = m_model->registerClass("IndependentSampler", "Independent sampler");
	m_rFilterNode = m_model->registerClass("BoxFilter", "Box filter");
	QRegExp resRegExp("^[1-9]\\d{0,4}x[1-9]\\d{0,4}$");
	ui->resolutionBox->setValidator(new QRegExpValidator(resRegExp, this));
	QPalette pal = ui->helpViewer->palette();
	pal.setColor(QPalette::Text, pal.color(QPalette::Foreground));
	pal.setColor(QPalette::Base, pal.color(QPalette::Window));
	ui->helpViewer->setPalette(pal);
}
	
void RenderSettingsDialog::setDocumentation(const QString &text) {
	m_currentDocumentation = text;
	bool hasErrors = false;
	QString comments;

	if (m_statusMessages.size() > 0) {
		comments = QString("<ul style=\"margin:0px;\">");
		ui->groupBox->setTitle(tr("Issues with the current configuration"));
		for (int i=0; i<m_statusMessages.size(); ++i) {
			const QString &message = m_statusMessages[i];
			bool isWarning = false, isError = false;

			if (message.contains("Warning"))
				isWarning = true;
			if (message.contains("Error"))
				isError = true;

			if (isError || isWarning)
				comments += QString("<li><em>%1</em></li>").arg(message);
			else
				comments += QString("<li>%1</li>").arg(message);
			hasErrors |= isError;
		}
		comments += "</ul>";
	} else {
		ui->groupBox->setTitle(tr("Documentation"));
	}

	ui->helpViewer->setHtml(comments + m_currentDocumentation);
   	ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(!hasErrors);
}

void RenderSettingsDialog::dataChanged() {
	QStringList statusMessages = validateConfiguration();
	if (statusMessages != m_statusMessages) {
		m_statusMessages = statusMessages;
		setDocumentation(m_currentDocumentation);
	}
}


void RenderSettingsDialog::cbHighlighted(int index) {
	QComboBox *comboBox = static_cast<QComboBox *>(sender());
	setDocumentation(comboBox->itemData(index).toList().at(1).toString());
}

void RenderSettingsDialog::chkBoxPressed() {
	QCheckBox *checkBox = static_cast<QCheckBox *>(sender());
	setDocumentation(checkBox->property("help").toString());
}

void RenderSettingsDialog::onTreeSelectionChange(const QItemSelection &selected, const QItemSelection &deselected) {
	QModelIndexList indexList = selected.indexes();
	if (indexList.size() > 0)
		setDocumentation(m_model->data(indexList[0], Qt::ToolTipRole).toString());
}

void RenderSettingsDialog::update() {
	int index = ui->integratorBox->currentIndex();

	int sampleCount = -1;
	if (sender() == ui->samplerBox) {
		Properties samplerProps;
		m_samplerNode->putProperties(samplerProps);
		if (samplerProps.hasProperty("sampleCount"))
			sampleCount = samplerProps.getInteger("sampleCount");
		else if (samplerProps.hasProperty("resolution"))
			sampleCount = std::pow(samplerProps.getInteger("resolution"), 2);
	}
	
	m_integratorNode = m_model->updateClass(m_integratorNode,
		ui->integratorBox->itemData(index).toList().at(0).toString(),
		ui->integratorBox->itemText(index));
	index = ui->samplerBox->currentIndex();
	m_samplerNode = m_model->updateClass(m_samplerNode,
		ui->samplerBox->itemData(index).toList().at(0).toString(),
		ui->samplerBox->itemText(index));
	index = ui->rFilterBox->currentIndex();
	m_rFilterNode = m_model->updateClass(m_rFilterNode,
		ui->rFilterBox->itemData(index).toList().at(0).toString(),
		ui->rFilterBox->itemText(index));

	if (ui->icBox->isChecked()) {
		m_icNode = m_model->updateClass(m_icNode,
			"IrradianceCacheIntegrator", tr("Irradiance Cache"));
	} else {
		m_icNode = m_model->updateClass(m_icNode, "", "");
	}

	if (ui->aiBox->isChecked()) {
		m_aiNode = m_model->updateClass(m_aiNode,
			"ErrorControl", tr("Adaptive Integration"));
	} else {
		m_aiNode = m_model->updateClass(m_aiNode, "", "");
	}

	if (sender() == ui->samplerBox && sampleCount != -1) {
		std::string samplerPlugin = getPluginName(ui->samplerBox);
		for (int i=0; i<m_samplerNode->childCount(); ++i) {
			TreeItem *treeItem = m_samplerNode->child(i);
			if (treeItem->getName() == "sampleCount") {
				treeItem->setValue(sampleCount);
			} else if (treeItem->getName() == "resolution") {
				treeItem->setValue((int) std::sqrt(sampleCount));
			}
		}
	}

	ui->treeView->expandAll();
	dataChanged();
}

bool RenderSettingsDialog::resolutionHasChanged() const {
	return ui->resolutionBox->currentText() != m_originalResolution;
}

void RenderSettingsDialog::refresh() {
	bool valid = true;
	int pos;

	QString resolutionString(ui->resolutionBox->currentText());
	valid &= ui->resolutionBox->validator()->validate(resolutionString,pos)
		== QValidator::Acceptable;

	ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(valid);
}

void RenderSettingsDialog::load(const SceneContext *ctx) {
	const Scene *scene = ctx->scene.get();
	const Film *film = scene->getFilm();
	const Properties &rFilterProps = film->getReconstructionFilter()->getProperties();
	const Properties &samplerProps = scene->getSampler()->getProperties();
	const Integrator *integrator = scene->getIntegrator();
	Properties integratorProps = integrator->getProperties();

	if (integratorProps.getPluginName() == "errctrl") {
		ui->aiBox->setChecked(true);
		m_model->setProperties(m_aiNode, integratorProps);
		integrator = integrator->getSubIntegrator();
		integratorProps = integrator->getProperties();
	}

	if (integratorProps.getPluginName() == "irrcache") {
		ui->icBox->setChecked(true);
		m_model->setProperties(m_icNode, integratorProps);
		integrator = integrator->getSubIntegrator();
		integratorProps = integrator->getProperties();
	}

	ui->resolutionBox->lineEdit()->setText(QString("%1x%2")
		.arg(film->getSize().x).arg(film->getSize().y));
	m_originalResolution = ui->resolutionBox->lineEdit()->text();

	setComboBox(ui->integratorBox, integratorProps.getPluginName());
	setComboBox(ui->rFilterBox, rFilterProps.getPluginName());
	setComboBox(ui->samplerBox, samplerProps.getPluginName());
	update();

	m_model->setProperties(m_rFilterNode, rFilterProps);
	m_model->setProperties(m_samplerNode, samplerProps);
	m_model->setProperties(m_integratorNode, integratorProps);
	ui->treeView->expandAll();
}

void RenderSettingsDialog::apply(SceneContext *ctx) {
	Scene *scene = new Scene(ctx->scene);
	const PinholeCamera *oldCamera = static_cast<const PinholeCamera *>(scene->getCamera());
	Properties filmProps = oldCamera->getFilm()->getProperties();
	ref<PluginManager> pluginMgr = PluginManager::getInstance();

	/* Configure the reconstruction filter */
	Properties rFilterProps(getPluginName(ui->rFilterBox));
	if (m_rFilterNode != NULL) 
		m_rFilterNode->putProperties(rFilterProps);
	ref<ReconstructionFilter> rFilter = static_cast<ReconstructionFilter *> 
		(pluginMgr->createObject(ReconstructionFilter::m_theClass, rFilterProps));

	/* Configure the sampler */
	Properties samplerProps(getPluginName(ui->samplerBox));
	if (m_samplerNode != NULL) 
		m_samplerNode->putProperties(samplerProps);
	ref<Sampler> sampler = static_cast<Sampler *> 
		(pluginMgr->createObject(Sampler::m_theClass, samplerProps));
	sampler->configure();
	
	/* Configure the integrator */
	Properties integratorProps(getPluginName(ui->integratorBox));
	if (m_integratorNode != NULL) 
		m_integratorNode->putProperties(integratorProps);
	ref<Integrator> integrator = static_cast<Integrator *> 
		(pluginMgr->createObject(Integrator::m_theClass, integratorProps));
	integrator->configure();

	if (ui->icBox->isChecked()) {
		Properties icProps("irrcache");
		if (m_icNode != NULL) 
			m_icNode->putProperties(icProps);
		ref<Integrator> ic = static_cast<Integrator *> 
			(pluginMgr->createObject(Integrator::m_theClass, icProps));
		ic->addChild("", integrator);
		ic->configure();
		integrator = ic;
	}

	if (ui->aiBox->isChecked()) {
		Properties aiProps("errctrl");
		if (m_aiNode != NULL) 
			m_aiNode->putProperties(aiProps);
		ref<Integrator> ai = static_cast<Integrator *> 
			(pluginMgr->createObject(Integrator::m_theClass, aiProps));
		ai->addChild("", integrator);
		ai->configure();
		integrator = ai;
	}

	/* Configure the film */
	QStringList resolution = ui->resolutionBox->currentText().split('x');
	SAssert(resolution.size() == 2);
	int width = resolution[0].toInt(), height = resolution[1].toInt();
	filmProps.setInteger("width", width, false);
	filmProps.setInteger("height", height, false);
	ref<Film> film = static_cast<Film *> (pluginMgr->createObject(Film::m_theClass, filmProps));
	film->addChild("", rFilter);
	film->configure();
	if (width != ctx->framebuffer->getWidth() ||
		height != ctx->framebuffer->getHeight()) {
		ctx->framebuffer = new Bitmap(width, height, 128);
		ctx->framebuffer->clear();
		ctx->mode = EPreview;
	}

	/* Configure the camera */
	Properties cameraProps = oldCamera->getProperties();
	ref<PinholeCamera> camera = static_cast<PinholeCamera *> 
		(pluginMgr->createObject(Camera::m_theClass, cameraProps));
	camera->addChild("", sampler);
	camera->addChild("", film);
	camera->setViewTransform(oldCamera->getViewTransform());
	camera->setFov(oldCamera->getFov());
	camera->configure();
	
	/* Update the scene with the newly constructed elements */
	scene->setCamera(camera);
	scene->setSampler(sampler);
	scene->setIntegrator(integrator);
	scene->configure();

	ctx->scene = scene;
}

RenderSettingsDialog::~RenderSettingsDialog() {
	delete ui;
}

void RenderSettingsDialog::changeEvent(QEvent *e) {
    QDialog::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
		ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

/* ====================== PropertyDelegate impl ====================== */

PropertyDelegate::PropertyDelegate(QObject *parent) : QStyledItemDelegate(parent) {
}

PropertyDelegate::~PropertyDelegate() {
}

QString	PropertyDelegate::displayText(const QVariant &value, const QLocale &locale) const {
	if (value.type() == QVariant::Bool)
		return value.toBool() ? tr("Yes") : tr("No");
	return QStyledItemDelegate::displayText(value, locale);
}

QWidget *PropertyDelegate::createEditor(QWidget *parent, const QStyleOptionViewItem &option,
		const QModelIndex &index) const {
	if (index.data().type() == QVariant::Bool) {
		QComboBox *cbox = new QComboBox(parent);
		/* Nicer boolean editor -- by default, Qt creates a True/False combo box */
		cbox->addItem(tr("No"));
		cbox->addItem(tr("Yes"));
		return cbox;
	}
	QWidget *widget = QStyledItemDelegate::createEditor(parent, option, index);
#if defined(__OSX__)
	/* Don't draw focus halos on OSX, they're really distracting */
	if (widget != NULL && widget->testAttribute(Qt::WA_MacShowFocusRect))
		widget->setAttribute(Qt::WA_MacShowFocusRect, false);
	if (index.data().type() != QVariant::Bool) {
		widget->setAttribute(Qt::WA_MacMiniSize, true);
		widget->setStyleSheet("font-size: 13pt;");
	}
#endif
	return widget;
}

void PropertyDelegate::setEditorData(QWidget *editor, const QModelIndex &index) const {
	if (index.data().type() == QVariant::Bool) {
		QComboBox *cbox = static_cast<QComboBox *>(editor);
		cbox->setCurrentIndex(index.data().toBool() ? 1 : 0);
		return;
	}

	QStyledItemDelegate::setEditorData(editor, index);
}

void PropertyDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
	const QModelIndex &index) const {
	if (index.data().type() == QVariant::Bool) {
		QComboBox *cbox = static_cast<QComboBox *>(editor);
		model->setData(index, QVariant(cbox->currentIndex() == 1), Qt::EditRole);
		return;
	}
	QStyledItemDelegate::setModelData(editor, model, index);
}

void PropertyDelegate::updateEditorGeometry(QWidget *editor,
	const QStyleOptionViewItem &option, const QModelIndex &index) const {
	if (index.data().type() == QVariant::Bool) {
		editor->setGeometry(option.rect);
		return;
	}

	QStyledItemDelegate::updateEditorGeometry(editor, option, index);
}

QStringList RenderSettingsDialog::validateConfiguration() const {
	/* Ad-hoc verification until we have something better 
	   (preferably specifiable by the plugins themselves) */
	QStringList messages;
	std::string integratorName = getPluginName(ui->integratorBox);
	std::string samplerName = getPluginName(ui->samplerBox);
	Properties integratorProps, samplerProps;
	m_integratorNode->putProperties(integratorProps);
	m_samplerNode->putProperties(samplerProps);

	if (samplerName != "independent") {
		if (integratorName == "ptracer" || integratorName == "mlt" || integratorName == "kelemen")
			messages << "Error: This integrator requires the independent sampler.";
		if (ui->aiBox->isChecked())
			messages << "Error: Adaptive integration requires the independent sampler.";
	}
	if (ui->icBox->isChecked()) {
		if (integratorName != "direct" && integratorName != "path" && integratorName != "volpath" && integratorName != "volpath_simple" && integratorName != "photonmapper")
			messages << "Error: Irradiance Caching requires a sampling-based integrator.";

	}
	if (ui->aiBox->isChecked()) {
		if (integratorName != "direct" && integratorName != "path" && integratorName != "volpath" && integratorName != "volpath_simple" && integratorName != "photonmapper")
			messages << "Error: Adaptive integration requires a sampling-based integrator.";
	}
	if (integratorName == "ppm") {
		if ((samplerName == "independent" || samplerName == "ldsampler") && samplerProps.hasProperty("sampleCount")) {
			if (samplerProps.getInteger("sampleCount") > 4)
				messages << "Warning: are you sure you need more than 4 samples/pixel for progressive photon mapping? This will be slow..";
		} else if (samplerName == "stratified" && samplerProps.hasProperty("resolution")) {
			if (samplerProps.getInteger("resolution") > 2)
				messages << "Warning: are you sure you need more than 4 samples/pixel for progressive photon mapping? This will be slow..";
		}
	}
	return messages;
}

