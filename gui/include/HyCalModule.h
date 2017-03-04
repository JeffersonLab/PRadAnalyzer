#ifndef QT_HYCAL_MODULE_H
#define QT_HYCAL_MODULE_H

#include <QGraphicsItem>
#include <QStyleOptionGraphicsItem>
#include <QFont>
#include "PRadHyCalModule.h"
#include "datastruct.h"


class PRadEventViewer;

class HyCalModule : public QGraphicsItem, public PRadHyCalModule
{
public:
    HyCalModule(PRadEventViewer* const p,
                const std::string &name,
                const PRadHyCalModule::Geometry &geo);
    virtual ~HyCalModule();

    void Initialize();
    void SetColor(const QColor &c) {color = c;};
    void SetColor(const double &val);
    void ShowPedestal();
    void ShowPedSigma();
    void ShowOccupancy();
    void ShowEnergy();
    void ShowEnergy(const double &energy);
    void ShowCustomValue() {SetColor(custom_value);};
    void SetCustomValue(double val) {custom_value = val;};
    const double &GetCustomValue() const {return custom_value;};
    QString GetReadID() const {return qname;};

    // overload
    QRectF boundingRect() const;
    void paint(QPainter *painter,
               const QStyleOptionGraphicsItem *option, QWidget *widget);
    void setSelected(bool selected);
    bool isSelected(){return m_selected;};

protected:
    void hoverEnterEvent(QGraphicsSceneHoverEvent *event);
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event);

private:
    PRadEventViewer *console;
    QString qname;

    bool m_hover;
    bool m_selected;
    QColor color;
    QFont font;
    QPainterPath shape;
    double custom_value;

#ifdef USE_CAEN_HV
public:
    void SetHVAddress(const ChannelAddress &set) {hv_addr = set;};
    const ChannelAddress &GetHVAddress() const {return hv_addr;};
private:
    ChannelAddress hv_addr;
#endif
};

#endif
