#ifndef PRAD_SPECTRUM_H
#define PRAD_SPECTRUM_H

#include <QGraphicsObject>
#include <QGradient>

class Spectrum : public QGraphicsObject
{
    Q_OBJECT

public:
    enum SpectrumType
    {
        Rainbow1,
        Rainbow2,
        GreyScale,
    };

    enum SpectrumScale
    {
        LogScale,
        LinearScale,
    };

    struct SettingData
    {
        SpectrumType type;
        SpectrumScale scale;
        double range_min;
        double range_max;

        SettingData()
        : type(Rainbow1), scale(LogScale), range_min(1.), range_max(100.)
        {};
    };

    Spectrum(double w, double h, double l1 = 440, double l2 = 640);
    void SetSpectrumRange(const double &x1, const double &x2);
    void SetSpectrumRangeMin(double &min);
    void SetSpectrumRangeMax(double &max);
    void SetSpectrumType(const SpectrumType &type);
    void SetSpectrumScale(const SpectrumScale &scale);
    SettingData &GetCurrentSetting() {return settings;};
    double GetRangeMin() {return settings.range_min;};
    double GetRangeMax() {return settings.range_max;};
    QColor GetColor(const double &val);
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    void paintTicks(QPainter *painter);
    QRectF boundingRect() const;

signals:
    void spectrumChanged();

private:
    void updateGradient();
    double scaling(const double &val);
    QColor scaleToColor(const double &scale);

private:
    SettingData settings;
    double width;
    double height;
    double wavelength1;
    double wavelength2;
    QLinearGradient gradient;
    QPainterPath shape;
};

#endif
