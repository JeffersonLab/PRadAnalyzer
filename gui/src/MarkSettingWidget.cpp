//============================================================================//
// A simple widget to set up the mark shape and color                         //
//                                                                            //
// Chao Peng                                                                  //
// 10/24/2016                                                                 //
//============================================================================//

#include "MarkSettingWidget.h"
#include <QCheckBox>
#include <QColorDialog>
#include <QHBoxLayout>
#include <QLabel>
#include <QComboBox>
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QDoubleSpinBox>

MarkSettingWidget::MarkSettingWidget(QWidget *parent)
: QWidget(parent), chosenColor(Qt::black)
{
    initialize();
}

MarkSettingWidget::MarkSettingWidget(const QStringList &m, QWidget *parent)
: QWidget(parent), savedColor(Qt::black), savedIndex(-1)
{
    initialize();
    SetMarkNames(m);
}

void MarkSettingWidget::AddMarkName(const QString &n)
{
    markCombo->addItem(n);
}

void MarkSettingWidget::initialize()
{
    checkBox = new QCheckBox("Show");
    checkBox->setChecked(true);

    matchCheckBox = new QCheckBox("Matched");
    matchCheckBox->setChecked(true);

    QLabel *markLabel = new QLabel(tr("Shape"));
    markLabel->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
    markCombo = new QComboBox;

    QLabel *markWidthLabel = new QLabel(tr("Width"));
    markWidth = new QSpinBox;
    markWidth->setRange(0, 20);

    QLabel *markSizeLabel = new QLabel(tr("Size"));
    markSizeLabel->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
    markSize = new QDoubleSpinBox;
    markSize->setRange(0., 20.0);
    markSize->setSingleStep(0.1);

    QLabel *colorLabel = new QLabel(tr("Color"));
    colorLabel->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
    colorButton = new QPushButton;
    colorButton->setMaximumWidth(30);
    colorButton->setMaximumHeight(25);

    QHBoxLayout *layout = new QHBoxLayout;
    layout->addWidget(checkBox);
    layout->addWidget(matchCheckBox);
    layout->addWidget(markWidthLabel);
    layout->addWidget(markWidth);
    layout->addWidget(markLabel);
    layout->addWidget(markCombo);
    layout->addWidget(markSizeLabel);
    layout->addWidget(markSize);
    layout->addWidget(colorLabel);
    layout->addWidget(colorButton);
    layout->setContentsMargins(0, 0, 0, 0);

    connect(colorButton, SIGNAL(clicked()), this, SLOT(colorPicker()));

    setLayout(layout);
}

bool MarkSettingWidget::IsChecked()
{
    return checkBox->isChecked();
}

bool MarkSettingWidget::IsMatchChecked()
{
    return matchCheckBox->isChecked();
}

void MarkSettingWidget::SetMarkNames(const QStringList &n)
{
    markCombo->clear();
    for(auto &name : n)
    {
        markCombo->addItem(name);
    }
}

void MarkSettingWidget::ClearMarkNames()
{
    markCombo->clear();
}

int MarkSettingWidget::GetCurrentMarkIndex()
{
    return markCombo->currentIndex();
}

QString MarkSettingWidget::GetCurrentMarkName()
{
    return markCombo->currentText();
}

int MarkSettingWidget::GetWidth()
{
    return markWidth->value();
}

QColor MarkSettingWidget::GetColor()
{
    return chosenColor;
}

double MarkSettingWidget::GetSize()
{
    return markSize->value();
}

void MarkSettingWidget::SetChecked(bool c)
{
    checkBox->setChecked(c);
}

void MarkSettingWidget::SetMatchChecked(bool c)
{
    matchCheckBox->setChecked(c);
}

void MarkSettingWidget::SetCurrentMarkIndex(int i)
{
    markCombo->setCurrentIndex(i);
}

void MarkSettingWidget::SetWidth(int w)
{
    markWidth->setValue(w);
}

void MarkSettingWidget::SetColor(const QColor &q)
{
    chosenColor = q;
    colorLabelChange();
}

void MarkSettingWidget::SetColor(const QString &c_text)
{
    chosenColor = QColor(c_text);
    colorLabelChange();
}

void MarkSettingWidget::SetSize(const double &s)
{
    markSize->setValue(s);
}

void MarkSettingWidget::SaveSettings()
{
    savedCheck = IsChecked();
    savedMatch = IsMatchChecked();
    savedIndex = GetCurrentMarkIndex();
    savedWidth = GetWidth();
    savedColor = GetColor();
    savedSize = GetSize();
}

void MarkSettingWidget::RestoreSettings()
{
    SetChecked(savedCheck);
    SetMatchChecked(savedMatch);
    SetCurrentMarkIndex(savedIndex);
    SetWidth(savedWidth);
    SetColor(savedColor);
    SetSize(savedSize);
}

void MarkSettingWidget::colorPicker()
{
    QColor color = QColorDialog::getColor(chosenColor, this, "Choose a color", QColorDialog::ShowAlphaChannel);

    if(!color.isValid())
        return;

    chosenColor = color;
    colorLabelChange();
}

void MarkSettingWidget::colorLabelChange()
{
    colorButton->setStyleSheet("background-color: rgba(" + QString::number(chosenColor.red())
                               + ", " + QString::number(chosenColor.green())
                               + ", " + QString::number(chosenColor.blue())
                               + ", " + QString::number(chosenColor.alpha()) + ");");
    colorButton->update();
}
