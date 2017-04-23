#ifndef PRAD_EVENT_VIEWER_H
#define PRAD_EVENT_VIEWER_H

#include <QMainWindow>
#include <QFileDialog>
#include <QFutureWatcher>
#include <vector>

#define HYCAL_SHIFT -50
#define CARTESIAN_TO_HYCALSCENE(x, y) x+HYCAL_SHIFT, -y

class HyCalScene;
class HyCalView;
class HyCalModule;
class Spectrum;
class SpectrumSettingPanel;
class HistCanvas;
class LogsBox;
class PRadDataHandler;
class PRadEPICSystem;
class PRadTaggerSystem;
class PRadHyCalSystem;
class PRadGEMSystem;
class PRadCoordSystem;
class PRadDetMatch;

#ifdef RECON_DISPLAY
class ReconSettingPanel;
#endif

#ifdef USE_ONLINE_MODE
class PRadETChannel;
class ETSettingPanel;
#endif

#ifdef USE_CAEN_HV
class PRadHVSystem;
#endif

QT_BEGIN_NAMESPACE
class QPushButton;
class QComboBox;
class QSpinBox;
class QSlider;
class QString;
class QLabel;
class QSplitter;
class QTreeWidget;
class QTreeWidgetItem;
class QTimer;
class QAction;
QT_END_NAMESPACE

enum HistType {
    EnergyTDCHist,
    ModuleHist,
    TaggerHist,
};

enum AnnoType {
    NoAnnotation,
    ShowID,
    ShowDAQ,
    ShowTDC,
};

enum ViewMode {
    EnergyView,
    OccupancyView,
    PedestalView,
    SigmaView,
    CustomView,
    HighVoltageView,
    VoltageSetView,
};

enum ViewerStatus {
    NO_INPUT,
    DATA_FILE,
    ONLINE_MODE,
};

class PRadEventViewer : public QMainWindow
{
    Q_OBJECT

public:
    PRadEventViewer();
    virtual ~PRadEventViewer();
    ViewMode GetViewMode() {return viewMode;}
    AnnoType GetAnnoType() {return annoType;}
    QColor GetColor(const double &val);
    void UpdateStatusBar(ViewerStatus mode);
    void UpdateStatusInfo();
    void UpdateHistCanvas();
    void SelectModule(HyCalModule* module);
    PRadDataHandler *GetHandler() {return handler;}

signals:
    void currentEventChanged(int evt);

public slots:
    void Refresh();

private slots:
    void handleEventChange(int event);
    void openDataFile();
    void initializeFromFile();
    void openCalibrationFile();
    void openGainFactorFile();
    void openCustomMap();
    void handleRootEvents();
    void saveHistToFile();
    void findPeak();
    void fitPedestal();
    void fitHistogram();
    void correctGainFactor();
    void takeSnapShot();
    void changeHistType(int index);
    void changeAnnoType(int index);
    void changeViewMode(int index);
    void changeSpectrumSetting();
    void changeCurrentEvent(int evt);
    void eraseBufferAction();
    void findEvent();
    void editCustomValueLabel(QTreeWidgetItem* item, int column);

private:
    void initView();
    void setupUI();
    void generateSpectrum();
    void generateHyCalModules();
    void generateScalerBoxes();
    void setTDCGroupBox();
    void eraseData();
    void createMainMenu();
    void createControlPanel();
    void createStatusBar();
    void createStatusWindow();
    void setupInfoWindow();
    void updateEventRange();
    void readEventFromFile(const QString &filepath);
    void readCustomValue(const QString &filepath);
    void onlineUpdate(const size_t &max_events);
    bool onlineSettings();
    QMenu *setupFileMenu();
    QMenu *setupCalibMenu();
    QMenu *setupToolMenu();
    QMenu *setupSettingMenu();
    QString getFileName(const QString &title,
                        const QString &dir,
                        const QStringList &filter,
                        const QString &suffix,
                        QFileDialog::AcceptMode mode = QFileDialog::AcceptOpen);
    QStringList getFileNames(const QString &title,
                             const QString &dir,
                             const QStringList &filter,
                             const QString &suffix,
                             QFileDialog::AcceptMode mode = QFileDialog::AcceptOpen,
                             QFileDialog::FileMode fmode = QFileDialog::ExistingFiles);

    PRadDataHandler *handler;
    PRadEPICSystem *epic_sys;
    PRadTaggerSystem *tagger_sys;
    PRadHyCalSystem *hycal_sys;
    PRadGEMSystem *gem_sys;
    HistType histType;
    AnnoType annoType;
    ViewMode viewMode;

    HyCalModule *selection;
    Spectrum *energySpectrum;
    //GEM *myGEM;
    HyCalScene *HyCal;
    HyCalView *view;
    HistCanvas *histCanvas;

    QString fileName;

    QSplitter *statusWindow;
    QSplitter *rightPanel;
    QSplitter *mainSplitter;

    QWidget *controlPanel;
    QSpinBox *eventSpin;
    QLabel *eventCntLabel;
    QComboBox *histTypeBox;
    QComboBox *annoTypeBox;
    QComboBox *viewModeBox;
    QPushButton *spectrumSettingButton;

    QTreeWidget *statusInfoWidget;
    QTreeWidgetItem *statusItem[6];

    QLabel *lStatusLabel;
    QLabel *rStatusLabel;

    QAction *openDataAction;

    QFileDialog *fileDialog;
    SpectrumSettingPanel *specSetting;
    LogsBox *logBox;

    QFuture<bool> future;
    QFutureWatcher<void> watcher;

#ifdef USE_ONLINE_MODE
public:
    void UpdateOnlineInfo();
private slots:
    void initOnlineMode();
    bool connectETClient();
    void startOnlineMode();
    void stopOnlineMode();
    void handleOnlineTimer();
private:
    void setupOnlineMode();
    QMenu *setupOnlineMenu();

    PRadETChannel *etChannel;
    QTimer *onlineTimer;
    ETSettingPanel *etSetting;
    QAction *onlineEnAction;
    QAction *onlineDisAction;
#endif

#ifdef USE_CAEN_HV
signals:
    void HVSystemInitialized();
private slots:
    void connectHVSystem();
    void initHVSystem();
    void disconnectHVSystem();
    void startHVMonitor();
    void saveHVSetting();
    void restoreHVSetting();
private:
    void setupHVSystem(const QString &list_file);
    QMenu *setupHVMenu();

    PRadHVSystem *hvSystem;
    QAction *hvEnableAction;
    QAction *hvDisableAction;
    QAction *hvSaveAction;
    QAction *hvRestoreAction;
#endif

#ifdef RECON_DISPLAY
private slots:
    void showReconEvent();
    void setupReconMethods();
    void enableReconstruct();
private:
    void setupReconDisplay();

    PRadCoordSystem *coordSystem;
    PRadDetMatch *detMatch;
    ReconSettingPanel *reconSetting;
    QSpinBox *clusterSpin;
#endif
};

#endif
