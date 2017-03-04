//============================================================================//
// Derived from QDialog, provide a setting panel to Change Cluster Display    //
//                                                                            //
// Chao Peng                                                                  //
// 10/24/2016                                                                 //
//============================================================================//

#include "ReconSettingPanel.h"
#include "PRadHyCalSystem.h"
#include "PRadGEMSystem.h"
#include "PRadCoordSystem.h"
#include "PRadDetMatch.h"
#include "MarkSettingWidget.h"

#include <QFileDialog>
#include <QVBoxLayout>
#include <QFormLayout>
#include <QGridLayout>
#include <QComboBox>
#include <QGroupBox>
#include <QLineEdit>
#include <QPushButton>
#include <QDialogButtonBox>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLabel>

#define COORD_ITEMS 6
#define MATCH_ITEMS 6

ReconSettingPanel::ReconSettingPanel(QWidget *parent)
: QDialog(parent), hycal(nullptr), gem(nullptr), coordSystem(nullptr), detMatch(nullptr)
{
    setWindowTitle(tr("Reconstruction Settings"));

    // set the pointer container size
    coordBox.resize(COORD_ITEMS, nullptr);
    matchConfLabel.resize(MATCH_ITEMS, nullptr);
    matchConfBox.resize(MATCH_ITEMS, nullptr);
    markSettings.resize(PRadDetector::Max_Dets, nullptr);

    QFormLayout *grid = new QFormLayout;
    grid->addRow(createMarkGroup());
    grid->addRow(createHyCalGroup(), createGEMGroup());
    grid->addRow(createCoordGroup());
    grid->addRow(createMatchGroup());
    grid->addRow(createStandardButtons());

    setLayout(grid);
}

QGroupBox *ReconSettingPanel::createMarkGroup()
{
    markGroup = new QGroupBox(tr("Cluster Mark Setting"));
    markGroup->setCheckable(true);
    markGroup->setChecked(true);
    QFormLayout *layout = new QFormLayout;

    // default settings
    int mark_index[] = {0, 4, 5};
    int mark_width[] = {2, 2, 2};
    double mark_size[] = {7.0, 5.0, 5.0};
    QColor mark_color[] = {Qt::black, Qt::red, Qt::magenta};

    for(size_t i = 0; i < markSettings.size(); ++i)
    {
        markSettings[i] = new MarkSettingWidget(HyCalScene::GetShapeList());
        markSettings[i]->SetCurrentMarkIndex(mark_index[i]);
        markSettings[i]->SetWidth(mark_width[i]);
        markSettings[i]->SetSize(mark_size[i]);
        markSettings[i]->SetColor(mark_color[i]);
        layout->addRow(tr(PRadDetector::getName(i)), markSettings[i]);
    }

    markGroup->setLayout(layout);

    return markGroup;
}

QGroupBox *ReconSettingPanel::createHyCalGroup()
{
    QGroupBox *hyCalGroup = new QGroupBox(tr("HyCal Cluster Setting"));

    // methods combo box
    hyCalMethods = new QComboBox;

    // configuration file
    QPushButton *hyCalLoadConfig = new QPushButton(tr("Reload"));
    QPushButton *hyCalFindPath = new QPushButton(tr("Open Configuration File"));
    hyCalConfigPath = new QLineEdit;

    connect(hyCalMethods, SIGNAL(currentIndexChanged(int)), this, SLOT(updateHyCalPath()));
    connect(hyCalLoadConfig, SIGNAL(clicked()), this, SLOT(loadHyCalConfig()));
    connect(hyCalFindPath, SIGNAL(clicked()), this, SLOT(openHyCalConfig()));

    QGridLayout *layout = new QGridLayout;
    layout->addWidget(hyCalMethods, 1, 0);
    layout->addWidget(hyCalFindPath, 1, 1, 1, 2);
    layout->addWidget(hyCalConfigPath, 2, 0, 1, 2);
    layout->addWidget(hyCalLoadConfig, 2, 2);

    hyCalGroup->setLayout(layout);

    return hyCalGroup;
}

QGroupBox *ReconSettingPanel::createGEMGroup()
{
    QGroupBox *gemGroup = new QGroupBox(tr("GEM Cluster Setting"));

    // the label name should be the same as configuration element name
    // it will be used to retrieve and set the configuration value
    gemMinLabel = new QLabel("Min Cluster Hits");
    gemMinHits = new QSpinBox;
    gemMinHits->setRange(0, 20);

    gemMaxLabel = new QLabel("Max Cluster Hits");
    gemMaxHits = new QSpinBox;
    gemMaxHits->setRange(2, 100);

    gemSplitLabel = new QLabel("Split Threshold");
    gemSplitThres = new QDoubleSpinBox;
    gemSplitThres->setRange(1., 1000.);
    gemSplitThres->setSingleStep(0.1);

    QFormLayout *layout = new QFormLayout;
    layout->addRow(gemMinLabel, gemMinHits);
    layout->addRow(gemMaxLabel, gemMaxHits);
    layout->addRow(gemSplitLabel, gemSplitThres);

    gemGroup->setLayout(layout);

    return gemGroup;
}

QGroupBox *ReconSettingPanel::createCoordGroup()
{
    QGroupBox *typeGroup = new QGroupBox(tr("Coordinates Setting"));

    coordType = new QComboBox;
    coordRun = new QComboBox;
    QPushButton *restoreCoord = new QPushButton("Restore");
    QPushButton *saveCoord = new QPushButton("Save");
    QPushButton *loadCoordFile = new QPushButton("Open Coord File");
    QPushButton *saveCoordFile = new QPushButton("Save Coord File");

    QGridLayout *layout = new QGridLayout;
    layout->addWidget(new QLabel(tr("Select Run")), 0, 0);
    layout->addWidget(coordRun, 0, 1);
    layout->addWidget(loadCoordFile, 0, 2);
    layout->addWidget(saveCoordFile, 0, 3);
    layout->addWidget(new QLabel(tr("Select Detector")), 1, 0);
    layout->addWidget(coordType, 1, 1);
    layout->addWidget(saveCoord, 1, 2);
    layout->addWidget(restoreCoord, 1, 3);
    layout->addWidget(new QLabel(tr("Axis")), 2, 0);
    layout->addWidget(new QLabel(tr("X")), 2, 1);
    layout->addWidget(new QLabel(tr("Y")), 2, 2);
    layout->addWidget(new QLabel(tr("Z")), 2, 3);
    layout->addWidget(new QLabel(tr("Origin (mm)")),3, 0);
    layout->addWidget(new QLabel(tr("Angle (mrad)")), 4, 0);

    int decimals[]     = {    4,     4,     1,     4,     4,     4};
    double step[]      = { 1e-4,  1e-4,  1e-1,  1e-4,  1e-4,  1e-4};
    double min_range[] = {  -20,   -20,     0,   -20,   -20,   -20};
    double max_range[] = {   20,    20, 10000,    20,    20,    20};

    for(size_t i = 0; i < coordBox.size(); ++i)
    {
        // setting coord spin box
        coordBox[i] = new QDoubleSpinBox;
        coordBox[i]->setDecimals(decimals[i]);
        coordBox[i]->setRange(min_range[i], max_range[i]);
        coordBox[i]->setSingleStep(step[i]);

        // add to layout
        int row = 3 + i/3;
        int col = 1 + i%3;
        layout->addWidget(coordBox[i], row, col);
    }

    typeGroup->setLayout(layout);

    connect(coordRun, SIGNAL(currentIndexChanged(int)), this, SLOT(selectCoordData(int)));
    connect(coordType, SIGNAL(currentIndexChanged(int)), this, SLOT(changeCoordType(int)));
    connect(restoreCoord, SIGNAL(clicked()), this, SLOT(restoreCoordData()));
    connect(saveCoord, SIGNAL(clicked()), this, SLOT(saveCoordData()));
    connect(loadCoordFile, SIGNAL(clicked()), this, SLOT(openCoordFile()));
    connect(saveCoordFile, SIGNAL(clicked()), this, SLOT(saveCoordFile()));

    return typeGroup;
}

QGroupBox *ReconSettingPanel::createMatchGroup()
{
    QGroupBox *matchGroup = new QGroupBox(tr("Coincidence Setting"));

    QGridLayout *layout = new QGridLayout;

    // the label name should be the same as configuration element name
    // it will be used to retrieve and set the configuration value
    QStringList matchDescript;
    matchDescript << "Lead Glass Resolution" << "Transition Resolution"
                  << "Crystal Resolution" << "Match Factor"
                  << "GEM Resolution" << "GEM Overlap Factor";

    for(size_t i = 0; i < matchConfBox.size(); ++i)
    {
        // label
        matchConfLabel[i] = new QLabel(matchDescript[i]);
        // conf spin box
        matchConfBox[i] = new QDoubleSpinBox;
        matchConfBox[i]->setRange(0., 50.);
        matchConfBox[i]->setSingleStep(0.01);
        // add to layout
        layout->addWidget(matchConfLabel[i], i/2, (i%2)*2);
        layout->addWidget(matchConfBox[i], i/2, (i%2)*2 + 1);
    }

    matchGroup->setLayout(layout);

    return matchGroup;
}

QDialogButtonBox *ReconSettingPanel::createStandardButtons()
{
    // Add standard buttons to layout
    QDialogButtonBox *buttonBox = new QDialogButtonBox(this);
    buttonBox->setStandardButtons(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);

    // Connect standard buttons
    connect(buttonBox->button(QDialogButtonBox::Ok), SIGNAL(clicked()),
                    this, SLOT(accept()));
    connect(buttonBox->button(QDialogButtonBox::Cancel), SIGNAL(clicked()),
                    this, SLOT(reject()));

    return buttonBox;
}

void ReconSettingPanel::ConnectHyCalSystem(PRadHyCalSystem *h)
{
    hycal = h;

    // set the methods combo box according to the data handler settings
    hyCalMethods->clear();

    if(h == nullptr)
        return;

    int index = -1; // means no value selected in combo box
    // name list
    auto methods = hycal->GetClusterMethodNames();
    // current method name
    auto current = hycal->GetClusterMethodName();

    for(size_t i = 0; i < methods.size(); ++i)
    {
        if(methods.at(i) == current)
            index = i;

        hyCalMethods->addItem(QString::fromStdString(methods.at(i)));
    }

    hyCalMethods->setCurrentIndex(index);
}

void ReconSettingPanel::ConnectGEMSystem(PRadGEMSystem *g)
{
    gem = g;

    if(g == nullptr)
        return;

    // set the config values from GEM clustering method
    PRadGEMCluster *gem_method = gem->GetClusterMethod();

    gemMinHits->setValue(gem_method->GetConfig<int>(gemMinLabel->text().toStdString()));
    gemMaxHits->setValue(gem_method->GetConfig<int>(gemMaxLabel->text().toStdString()));
    gemSplitThres->setValue(gem_method->GetConfig<double>(gemSplitLabel->text().toStdString()));
}

void ReconSettingPanel::ConnectCoordSystem(PRadCoordSystem *c)
{
    coordSystem = c;

    coordRun->clear();
    coordType->clear();

    if(c == nullptr)
        return;

    for(auto &it : c->GetCoordsData())
    {
        coordRun->addItem(QString::number(it.first));
    }

    det_coords = c->GetCurrentCoords();

    for(auto &coord : det_coords)
    {
        coordType->addItem(PRadDetector::getName(coord.det_enum));
    }

    coordType->setCurrentIndex(0);
}

void ReconSettingPanel::ConnectMatchSystem(PRadDetMatch *m)
{
    detMatch = m;

    for(size_t i = 0; i < matchConfBox.size(); ++i)
    {
        float value = 0.;
        if(detMatch != nullptr)
            value = detMatch->GetConfig<float>(matchConfLabel[i]->text().toStdString());

        matchConfBox[i]->setValue(value);
    }
}

void ReconSettingPanel::updateHyCalPath()
{
    std::string name = hyCalMethods->currentText().toStdString();
    PRadHyCalCluster *method = hycal->GetClusterMethod(name);

    if(method)
        hyCalConfigPath->setText(QString::fromStdString(method->GetConfigPath()));
    else
        hyCalConfigPath->setText("");
}

// open the configuration file for selected method
void ReconSettingPanel::openHyCalConfig()
{
    QString path = QFileDialog::getOpenFileName(this,
                                                "HyCal Cluster Configuration File",
                                                "config/",
                                                "text config file (*.conf *.txt)");

    if(path.isEmpty())
        return;

    hyCalConfigPath->setText(path);
    loadHyCalConfig();
}

// load the configuration file for selected method
void ReconSettingPanel::loadHyCalConfig()
{
    std::string method_name = hyCalMethods->currentText().toStdString();
    std::string config_path = hyCalConfigPath->text().toStdString();

    PRadHyCalCluster *method = hycal->GetClusterMethod(method_name);
    if(method)
        method->Configure(config_path);
}

void ReconSettingPanel::changeCoordType(int t)
{
    if(coordSystem == nullptr)
        return;

    if((size_t)t > coordBox.size()) { // -1 is the default nohting to select value
        for(auto &box : coordBox)
            box->setValue(0);
        return;
    }

    auto &coord = det_coords[t];

    for(size_t i = 0; i < coordBox.size(); ++i)
    {
        coordBox[i]->setValue(coord.get_dim_coord(i));
    }

}

void ReconSettingPanel::selectCoordData(int r)
{
    if(r < 0 || coordSystem == nullptr)
        return;

    coordSystem->ChooseCoord(r);
    det_coords = coordSystem->GetCurrentCoords();
    changeCoordType(coordType->currentIndex());
}

void ReconSettingPanel::saveCoordData()
{
    size_t idx = (size_t)coordType->currentIndex();
    if(idx >= det_coords.size())
        return;

    auto &coord = det_coords[idx];

    for(size_t i = 0; i < coordBox.size(); ++i)
    {
        coord.set_dim_coord(i, coordBox[i]->value());
    }
}

void ReconSettingPanel::saveCoordFile()
{
    QString path = QFileDialog::getSaveFileName(this,
                                                "Coordinates File",
                                                "config/",
                                                "text data file (*.dat *.txt)");
    if(path.isEmpty())
        return;

    saveCoordData();
    coordSystem->SetCurrentCoord(det_coords);
    coordSystem->SaveCoordData(path.toStdString());
}

void ReconSettingPanel::openCoordFile()
{
    QString path = QFileDialog::getOpenFileName(this,
                                                "Coordinates File",
                                                "config/",
                                                "text data file (*.dat *.txt)");
    if(path.isEmpty())
        return;

    coordSystem->LoadCoordData(path.toStdString());

    // re-connect to apply changes
    ConnectCoordSystem(coordSystem);
}

void ReconSettingPanel::restoreCoordData()
{
    if(coordSystem == nullptr)
        return;

    det_coords = coordSystem->GetCurrentCoords();
    changeCoordType(coordType->currentIndex());
}

// reconnects all the objects to read their settings
void ReconSettingPanel::SyncSettings()
{
    ConnectHyCalSystem(hycal);
    ConnectGEMSystem(gem);
    ConnectCoordSystem(coordSystem);
    ConnectMatchSystem(detMatch);
}

// save the settings that cannot by sync
void ReconSettingPanel::SaveSettings()
{
    for(auto &mark : markSettings)
        mark->SaveSettings();
}

// restore the saved settings
void ReconSettingPanel::RestoreSettings()
{
    for(auto &mark : markSettings)
        mark->RestoreSettings();
}

// apply all the changes
void ReconSettingPanel::ApplyChanges()
{
    // set hycal
    if(hycal) {
        // change HyCal clustering method
        std::string method = hyCalMethods->currentText().toStdString();
        hycal->SetClusterMethod(method);
    }

    // set gem
    if(gem) {
        // set corresponding gem cluster configuration values
        PRadGEMCluster *gem_method = gem->GetClusterMethod();
        gem_method->SetConfigValue(gemMinLabel->text().toStdString(), gemMinHits->value());
        gem_method->SetConfigValue(gemMaxLabel->text().toStdString(), gemMaxHits->value());
        gem_method->SetConfigValue(gemSplitLabel->text().toStdString(), gemSplitThres->value());
        // reaload the configuration
        gem_method->Configure();
    }

    // set coordinate system
    if(coordSystem) {
        saveCoordData();
        coordSystem->SetCurrentCoord(det_coords);
    }

    // set detector matching system
    if(detMatch) {
        for(size_t i = 0; i < matchConfBox.size(); ++i)
        {
            detMatch->SetConfigValue(matchConfLabel[i]->text().toStdString(), matchConfBox[i]->value());
        }
    }
}

bool ReconSettingPanel::IsEnabled()
const
{
    return markGroup->isChecked();
}

bool ReconSettingPanel::ShowDetector(int det)
const
{
    if(det < 0 || det >= PRadDetector::Max_Dets)
        return false;

    return markSettings[det]->IsChecked();
}

bool ReconSettingPanel::ShowMatchedDetector(int det)
const
{
    if(det < 0 || det >= PRadDetector::Max_Dets)
        return false;

    return markSettings[det]->IsMatchChecked();
}

HyCalScene::MarkAttributes ReconSettingPanel::GetMarkAttributes(int det)
const
{
    return HyCalScene::MarkAttributes(GetMarkIndex(det),
                                      GetMarkWidth(det),
                                      GetMarkColor(det),
                                      GetMarkSize(det));
}

int ReconSettingPanel::GetMarkIndex(int det)
const
{
    if(det < 0 || det >= PRadDetector::Max_Dets)
        return -1;

    return markSettings[det]->GetCurrentMarkIndex();
}

QString ReconSettingPanel::GetMarkName(int det)
const
{
    if(det < 0 || det >= PRadDetector::Max_Dets)
        return "";

    return markSettings[det]->GetCurrentMarkName();
}

int ReconSettingPanel::GetMarkWidth(int det)
const
{
    if(det < 0 || det >= PRadDetector::Max_Dets)
        return 1;

    return markSettings[det]->GetWidth();
}

QColor ReconSettingPanel::GetMarkColor(int det)
const
{
    if(det < 0 || det >= PRadDetector::Max_Dets)
        return Qt::black;

    return markSettings[det]->GetColor();
}

double ReconSettingPanel::GetMarkSize(int det)
const
{
    if(det < 0 || det >= PRadDetector::Max_Dets)
        return 0;

    return markSettings[det]->GetSize();
}
