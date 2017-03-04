//============================================================================//
// Derived from QDialog, provide a setting panel to change ET settings in GUI //
//                                                                            //
// Chao Peng                                                                  //
// 02/17/2016                                                                 //
//============================================================================//

#include "ETSettingPanel.h"
#include <QDialogButtonBox>
#include <QPushButton>
#include <QLabel>
#include <QFormLayout>
#include <QSpinBox>
#include <QLineEdit>

ETSettingPanel::ETSettingPanel(QWidget *parent)
: QDialog(parent)
{
    this->setWindowTitle("ET Channel Settings");
    QFormLayout *dialogLayout = new QFormLayout(this);

    QLabel *warnLabel = new QLabel("Do NOT change default settings unless you "
                                   "know what you are doing.");

    // Add settings for online mode
    QLabel *ipLabel = new QLabel("ET Host/Port");
    ipEdit = new QLineEdit(this);
    ipEdit->setText("clondaq6.jlab.org");

    portEdit = new QSpinBox(this);
    portEdit->setRange(1, 65535);
    portEdit->setValue(11111);

    QHBoxLayout *hostLayout = new QHBoxLayout();
    hostLayout->addWidget(ipEdit);
    hostLayout->addWidget(portEdit);

    QLabel *fileLabel = new QLabel("ET File");
    fileEdit = new QLineEdit(this);
    fileEdit->setText("/tmp/et_sys_clasprad");

    QLabel *stationLabel = new QLabel("Station Name");
    stationEdit = new QLineEdit(this);
    stationEdit->setText("online monitor");

    dialogLayout->addRow(warnLabel);
    dialogLayout->addRow(ipLabel, hostLayout);
    dialogLayout->addRow(fileLabel, fileEdit);
    dialogLayout->addRow(stationLabel, stationEdit);

    // Add standard buttons to layout
    QDialogButtonBox *buttonBox = new QDialogButtonBox(this);
    buttonBox->setStandardButtons(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    dialogLayout->addRow(buttonBox);

    // Connect standard buttons
    connect(buttonBox->button(QDialogButtonBox::Ok), SIGNAL(clicked()),
                    this, SLOT(accept()));
    connect(buttonBox->button(QDialogButtonBox::Cancel), SIGNAL(clicked()),
                    this, SLOT(reject()));
}

QString ETSettingPanel::GetETHost()
{
    return ipEdit->text();
}

int ETSettingPanel::GetETPort()
{
    return portEdit->value();
}

QString ETSettingPanel::GetETFilePath()
{
    return fileEdit->text();
}

QString ETSettingPanel::GetStationName()
{
    return stationEdit->text();
}
