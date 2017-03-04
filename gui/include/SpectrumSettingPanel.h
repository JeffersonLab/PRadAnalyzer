#ifndef SPECTRUM_SETTING_PANEL_H
#define SPECTRUM_SETTING_PANEL_H

#include <QDialog>


class QDoubleSpinBox;
class QSlider;
class QGroupBox;
class QRadioButton;
class Spectrum;

class SpectrumSettingPanel : public QDialog
{
    Q_OBJECT

public:
    SpectrumSettingPanel(QWidget *parent = 0);
    ~SpectrumSettingPanel() {};
    void ConnectSpectrum(Spectrum *s);
    void ChoosePreSetting(int val);
    void SetSpectrumRange(double min, double max);
    void SetLinearScale();
    void SetLogScale();

private slots:
    void changeScale();
    void changeType();
    void changeRangeMax(double value);
    void changeRangeMin(double value);
    void changePreSetting();

private:
    Spectrum *spectrum;
    bool preset;

    QGroupBox *createScaleGroup();
    QGroupBox *createTypeGroup();
    QGroupBox *createRangeGroup();
    QGroupBox *createPreSetGroup();

    QRadioButton *linearScale;
    QRadioButton *logScale;

    QRadioButton *rainbow1;
    QRadioButton *rainbow2;
    QRadioButton *greyscale;

    QRadioButton *energyView;
    QRadioButton *occupancyView;
    QRadioButton *voltageView;
    QRadioButton *pedestalView;
    QRadioButton *sigmaView;
    QRadioButton *customView;

    QDoubleSpinBox *minSpin;
    QDoubleSpinBox *maxSpin;
};

#endif
