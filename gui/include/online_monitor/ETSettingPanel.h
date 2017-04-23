#ifndef ET_SETTING_PANEL_H
#define ET_SETTING_PANEL_H

#include <QDialog>

class QLineEdit;
class QSpinBox;

class ETSettingPanel : public QDialog {
    Q_OBJECT
public:
    ETSettingPanel(QWidget *parent = 0);
    ~ETSettingPanel() {}
    QString GetETHost();
    QString GetETFilePath();
    QString GetStationName();
    int GetETPort();

private:
    QLineEdit *ipEdit;
    QSpinBox *portEdit;
    QLineEdit *fileEdit;
    QLineEdit *stationEdit;

};

#endif
