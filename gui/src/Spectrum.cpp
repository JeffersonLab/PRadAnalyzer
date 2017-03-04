//============================================================================//
// Derived from QGraphicsItem, a class of spectrum                            //
//                                                                            //
// Chao Peng                                                                  //
// 02/27/2016                                                                 //
//============================================================================//

#include <QPainter>
#include <QVector>
#include <algorithm>
#include <cmath>
#include "Spectrum.h"

Spectrum::Spectrum(double w, double h, double l1, double l2)
: width(w), height(h), wavelength1(l1), wavelength2(l2)
{
    shape.addRect(-width/2., -height/2., width, height);
    gradient = QLinearGradient(QPoint(width/2., height/2.),
                               QPoint(width/2., -height/2.));
    updateGradient();
}

void Spectrum::updateGradient()
{
    double scale;
    for(int i = 0; i < 100; ++ i)
    {
        scale = (double)i/100;
        gradient.setColorAt(scale, scaleToColor(scale));
    }
}

void Spectrum::SetSpectrumType(const SpectrumType &type)
{
    if(settings.type == type)
        return;

    settings.type = type;
    updateGradient();
    emit spectrumChanged();
}

void Spectrum::SetSpectrumScale(const SpectrumScale &scale)
{
    if(settings.scale == scale)
        return;

    settings.scale = scale;
    emit spectrumChanged();
}

void Spectrum::SetSpectrumRange(const double &x1, const double &x2)
{
    if(settings.range_min == x1 && settings.range_max == x2)
        return;

    settings.range_min = x1;
    settings.range_max = x2;
    emit spectrumChanged();
}

void Spectrum::SetSpectrumRangeMin(double &min)
{
    if(settings.range_min == min)
        return;

    settings.range_min = min;
    emit spectrumChanged();
}

void Spectrum::SetSpectrumRangeMax(double &max)
{
    if(settings.range_max == max)
        return;

    settings.range_max = max;
    emit spectrumChanged();
}

// get color in log scale
QColor Spectrum::GetColor(const double &val)
{
    if(val < settings.range_min || (settings.range_min == settings.range_max))
        return Qt::white;
    else
        return scaleToColor(scaling(val));
}

double Spectrum::scaling(const double &val)
{
    double scale, rmin, rmax;
    if(settings.scale == LogScale) {
        rmin = log10(std::max(0.1, settings.range_min));
        rmax = log10(std::max(1.0, settings.range_max));
        scale = (log10(val) - rmin)/(rmax - rmin);
    } else {
        rmin = (double)settings.range_min;
        rmax = (double)settings.range_max;
        scale = (val - rmin)/(rmax - rmin);
    }

    if(scale < 0)
        scale = 0;
    else if(scale > 1)
        scale = 1;

    return scale;
}

// calculte the color
QColor Spectrum::scaleToColor(const double &scale)
{
    // Based on: http://www.efg2.com/Lab/ScienceAndEngineering/Spectra.htm
    // The foregoing is based on: http://www.midnightkite.com/color.html
    //double lambda = 1.0/(1.0/wavelength1 - scale*(1.0/wavelength1 - 1.0/wavelength2));
    double lambda = wavelength1 + scale*(wavelength2 - wavelength1);

    QColor rainbowColor;

    switch(settings.type) {
        case Rainbow1: // raibown type 1
        {
            struct Color {
                qreal c[3];
                QColor toColor(qreal f) const {
                    qreal const gamma = 0.8;
                    int ci[3];
                    for (int i = 0; i < 3; ++i) {
                       ci[i] = c[i] == 0.0 ? 0.0 : qRound(255 * pow(c[i] * f, gamma));
                    }
                    return QColor(ci[0], ci[1], ci[2]);
                }
            } color;

            qreal factor = 0.0;
            color.c[0] = color.c[1] = color.c[2] = 0.0;
            static qreal thresholds[] = { 380, 440, 490, 510, 580, 645, 780 };
            for (unsigned int i = 0; i < sizeof(thresholds)/sizeof(thresholds[0]); ++ i) {
                qreal t1 = thresholds[i], t2 = thresholds[i+1];
                if (lambda < t1 || lambda >= t2) continue;
                if (i&1) std::swap(t1, t2);
                color.c[i%3] = (i < 5) ? (lambda - t2) / (t1-t2) : 0.0;;
                color.c[2-i/2] = 1.0;
                factor = 1.0;
                break;
            }

            // Let the intensity fall off near the vision limits
            if (lambda >= 380 && lambda < 420) {
                factor = 0.3 + 0.7*(lambda-380) / (420 - 380);
            }
            else if (lambda >= 700 && lambda < 780) {
                factor = 0.3 + 0.7*(780 - lambda) / (780 - 700);
            }
            rainbowColor = color.toColor(factor);
        }
        break;

        case Rainbow2: // rainbow type 2
        {
            double l = lambda;
            double t,  r=0.0, g=0.0, b=0.0;
                 if ((l>=400.0)&&(l<410.0)) { t=(l-400.0)/(410.0-400.0); r=    +(0.33*t)-(0.20*t*t); }
            else if ((l>=410.0)&&(l<475.0)) { t=(l-410.0)/(475.0-410.0); r=0.14         -(0.13*t*t); }
            else if ((l>=545.0)&&(l<595.0)) { t=(l-545.0)/(595.0-545.0); r=    +(1.98*t)-(     t*t); }
            else if ((l>=595.0)&&(l<650.0)) { t=(l-595.0)/(650.0-595.0); r=0.98+(0.06*t)-(0.40*t*t); }
            else if ((l>=650.0)&&(l<700.0)) { t=(l-650.0)/(700.0-650.0); r=0.65-(0.84*t)+(0.20*t*t); }
                 if ((l>=415.0)&&(l<475.0)) { t=(l-415.0)/(475.0-415.0); g=             +(0.80*t*t); }
            else if ((l>=475.0)&&(l<590.0)) { t=(l-475.0)/(590.0-475.0); g=0.8 +(0.76*t)-(0.80*t*t); }
            else if ((l>=585.0)&&(l<639.0)) { t=(l-585.0)/(639.0-585.0); g=0.84-(0.84*t)           ; }
                 if ((l>=400.0)&&(l<475.0)) { t=(l-400.0)/(475.0-400.0); b=    +(2.20*t)-(1.50*t*t); }
            else if ((l>=475.0)&&(l<560.0)) { t=(l-475.0)/(560.0-475.0); b=0.7 -(     t)+(0.30*t*t); }

            rainbowColor = QColor(r*255, g*255, b*255);
        }
        break;

        case GreyScale: // gray scale
        {
            rainbowColor = QColor((1.0 - scale)*255,(1.0 - scale)*255,(1.0 - scale)*255);
        }
        break;
    }

    return rainbowColor;
}

// paint the spectrum
void Spectrum::paint(QPainter *painter,
                     const QStyleOptionGraphicsItem * /*option*/, QWidget * /* widget */)
{
    painter->fillPath(shape, gradient);
    paintTicks(painter);
}

void Spectrum::paintTicks(QPainter *painter)
{
    QVector<double> ticks;
    double scale, left;
    painter->setPen(Qt::black);
    painter->setFont(QFont("Courier",14, QFont::Bold));

    // calculate ticks
    switch(settings.scale)
    {
    case LogScale:
        ticks.push_back(std::max(0.1, settings.range_min));
        for(int i = 0; ; ++i)
        {
            int tick = pow(10, i);
            if(tick >= settings.range_max)
                break;
            if(tick > settings.range_min)
                ticks.push_back(tick);
        }
        ticks.push_back(std::max(1.0, settings.range_max));

        // if ticks are too less
        if(ticks.size() <= 3) {
            int ori_size = ticks.size();
            for(int i = 0; i < ori_size; ++i)
            {
                for(int j = 2; j < 10; ++j)
                {
                    if(j*ticks.at(i) >= settings.range_max)
                        break;
                    ticks.push_back(j*ticks.at(i));
                }
            }
        }
        else if(ticks.size() <= 5) {
            int ori_size = ticks.size();
            for(int i = 0; i < ori_size; ++i)
            {
               if(5*ticks.at(i) >= settings.range_max)
                   break;
               ticks.push_back(5*ticks.at(i));
            }
        }
        break;
    case LinearScale:
        ticks.push_back(settings.range_min);
        double step = (double)(settings.range_max - settings.range_min)/10.;
        for(int i = 0; i <= 10; ++i)
        {
            double tick = settings.range_min + step * i;
            if(tick != ticks.last())
                ticks.push_back(tick);
        }
        break;
    }

    for(int i = 0; i < ticks.size(); ++i)
    {
        scale = scaling(ticks[i]);
        if(scale > 0 && scale < 1)
            left = 0.;
        else
            left = -width/2.;
        painter->drawLine(left, (0.5 - scale)*height, width/2. + 50, (0.5 - scale)*height);
        QRect textBox = QRect(width/2 + 2., (0.5 - scale)*height - 22, 80, 20);
        painter->drawText(textBox, QString::number(ticks[i]), QTextOption(Qt::AlignBottom | Qt::AlignLeft));
    }
}

// define the bound of this item
QRectF Spectrum::boundingRect() const
{
    return QRectF(-width/2., -height/2., width, height);
}
