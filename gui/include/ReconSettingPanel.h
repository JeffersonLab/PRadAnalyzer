#ifndef RECON_SETTING_PANEL_H
#define RECON_SETTING_PANEL_H

#include <string>
#include <vector>
#include <QDialog>
#include "HyCalScene.h"
#include "PRadCoordSystem.h"
#include "PRadDetector.h"

class PRadHyCalSystem;
class PRadGEMSystem;
class PRadDetMatch;
class MarkSettingWidget;

class QLabel;
class QSpinBox;
class QDoubleSpinBox;
class QComboBox;
class QGroupBox;
class QLineEdit;
class QDialogButtonBox;

class ReconSettingPanel : public QDialog
{
    Q_OBJECT

public:
    ReconSettingPanel(QWidget *parent = 0);
    ~ReconSettingPanel() {};

    // connect objects
    void ConnectHyCalSystem(PRadHyCalSystem *h);
    void ConnectGEMSystem(PRadGEMSystem *g);
    void ConnectCoordSystem(PRadCoordSystem *c);
    void ConnectMatchSystem(PRadDetMatch *m);

    // change related
    void SyncSettings();
    void SaveSettings();
    void RestoreSettings();
    void ApplyChanges();

    // get data
    bool IsEnabled() const;
    bool ShowDetector(int det) const;
    bool ShowMatchedDetector(int det) const;
    HyCalScene::MarkAttributes GetMarkAttributes(int det) const;
    int GetMarkIndex(int det) const;
    QString GetMarkName(int det) const;
    int GetMarkWidth(int det) const;
    QColor GetMarkColor(int det) const;
    double GetMarkSize(int det) const;

private:
    QGroupBox *createMarkGroup();
    QGroupBox *createHyCalGroup();
    QGroupBox *createGEMGroup();
    QGroupBox *createCoordGroup();
    QGroupBox *createMatchGroup();
    QDialogButtonBox *createStandardButtons();

private slots:
    void updateHyCalPath();
    void loadHyCalConfig();
    void openHyCalConfig();
    void selectCoordData(int r);
    void changeCoordType(int t);
    void saveCoordData();
    void saveCoordFile();
    void openCoordFile();
    void restoreCoordData();

private:
    PRadHyCalSystem *hycal;
    PRadGEMSystem *gem;
    PRadCoordSystem *coordSystem;
    PRadDetMatch *detMatch;

    QGroupBox *markGroup;
    QComboBox *hyCalMethods;
    QLineEdit *hyCalConfigPath;

    QLabel *gemMinLabel;
    QSpinBox *gemMinHits;
    QLabel *gemMaxLabel;
    QSpinBox *gemMaxHits;
    QLabel *gemSplitLabel;
    QDoubleSpinBox *gemSplitThres;

    QComboBox *coordRun;
    QComboBox *coordType;
    std::vector<QDoubleSpinBox*> coordBox; // X, Y, Z, thetaX, thetaY, thetaZ

    std::vector<QLabel*> matchConfLabel;
    std::vector<QDoubleSpinBox*> matchConfBox; // resolution and matching factors

    std::vector<MarkSettingWidget*> markSettings;

private:
    std::vector<PRadCoordSystem::DetCoord> det_coords;
};

#endif
