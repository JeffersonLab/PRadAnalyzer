//============================================================================//
// HyCal Module class, the basic element of HyCal                             //
//                                                                            //
// Chao Peng                                                                  //
// 02/27/2016                                                                 //
//============================================================================//

#include <cmath>
#include <string>
#include "HyCalModule.h"
#include "PRadEventViewer.h"
#include "PRadADCChannel.h"
#if QT_VERSION >= 0x500000
#include <QtWidgets>
#else
#include <QtGui>
#endif

HyCalModule::HyCalModule(PRadEventViewer* const p,
                         const std::string &n,
                         const PRadHyCalModule::Geometry &geo)
: PRadHyCalModule(n, geo),
  console(p), qname(QString::fromStdString(n)), m_hover(false),
  m_selected(false), color(Qt::white), font(QFont("times",10)), custom_value(0.)
{
    // initialize the item
    Initialize();

    // detect if mouse is hovering on this item
    setAcceptHoverEvents(true);
}

HyCalModule::~HyCalModule()
{
    // place holder
}

// initialize the module
void HyCalModule::Initialize()
{
    // flip Y axis to transform Cartesian coords to QtWidget Coords
    // add HYCAL_SHIFT to make the HyCal display in center
    QGraphicsItem::setPos(CARTESIAN_TO_HYCALSCENE(geometry.x, geometry.y));

    if(geometry.type == PbWO4)
        font.setPixelSize(9);

    // graphical shape
    double size_x = geometry.size_x;
    double size_y = geometry.size_y;
    shape.addRect(-size_x/2., -size_y/2., size_x, size_y);
}

// define the bound of this item
QRectF HyCalModule::boundingRect()
const
{
    double size_x = geometry.size_x;
    double size_y = geometry.size_y;
    return QRectF(-size_x/2., -size_y/2., size_x, size_y);
}

// how to paint this item
void HyCalModule::paint(QPainter *painter,
                        const QStyleOptionGraphicsItem *option,
                        QWidget * /* widget */)
{
    QColor fillColor = color;
    QColor fontColor = Qt::black;

    // mouse is hovering around, turns darker
    if(m_hover == true) {
        fillColor = fillColor.darker(125);
        fontColor = Qt::white;
    }

    // if it is zoomed in a lot, put some gradient on the color
    if(option->levelOfDetail < 4.0) {
        painter->fillPath(shape, fillColor);
    } else {
        QLinearGradient gradient(QPoint(-20, -20), QPoint(+20, +20));
        int coeff = 105 + int(std::log(option->levelOfDetail - 4.0));
        gradient.setColorAt(0.0, fillColor.lighter(coeff));
        gradient.setColorAt(1.0, fillColor.darker(coeff));
        painter->fillPath(shape, gradient);
    }

    // draw frame of this item
    painter->drawPath(shape);

    // is it selected? turns purple if so
    if(m_selected) {
        painter->fillPath(shape,QColor(155,0,155,200));
        fontColor = Qt::white;
    }

    // get the current display settings from console
    if(console) {
        switch(console->GetAnnoType())
        {
        case NoAnnotation: // show nothing
            break;
        case ShowID: // show module id
            painter->setFont(font);
            painter->setPen(fontColor);
            if(geometry.type == PbWO4) // ignore W for PbWO4 module due to have a better view
                painter->drawText(boundingRect(), qname.mid(1), QTextOption(Qt::AlignCenter));
            else
                painter->drawText(boundingRect(), qname, QTextOption(Qt::AlignCenter));
            break;
        case ShowDAQ: // show daq configuration
            if(daq_ch) {
                unsigned int crate = daq_ch->GetAddress().crate;
                unsigned int slot = daq_ch->GetAddress().slot;
                unsigned int channel = daq_ch->GetAddress().channel;
                painter->setFont(font);
                painter->setPen(fontColor);
                painter->drawText(boundingRect(),
                              QString::number(crate) + "," + QString::number(slot),
                              QTextOption(Qt::AlignTop | Qt::AlignHCenter));
                painter->drawText(boundingRect(),
                              "Ch" + QString::number(channel+1),
                              QTextOption(Qt::AlignBottom | Qt::AlignHCenter));
            } else {
                painter->drawText(boundingRect(), "N/A", QTextOption(Qt::AlignCenter));
            }
            break;
        case ShowTDC: // this will be handled in HyCalScene
            break;
        }
    }
}

// hover bool
void HyCalModule::hoverEnterEvent(QGraphicsSceneHoverEvent *)
{
    m_hover = true;
    update();
}

void HyCalModule::hoverLeaveEvent(QGraphicsSceneHoverEvent *)
{
    m_hover = false;
    update();
}

void HyCalModule::setSelected(bool selected)
{
    m_selected = selected;
    update();
    if(selected && console)
        console->SelectModule(this);
}

// Get color from the spectrum
void HyCalModule::SetColor(const double &val)
{
    if(console)
        color = console->GetColor(val);
}

void HyCalModule::ShowPedestal()
{
    if(daq_ch)
        SetColor(daq_ch->GetPedestal().mean);
}

void HyCalModule::ShowPedSigma()
{
    if(daq_ch)
        SetColor(daq_ch->GetPedestal().sigma);
}

void HyCalModule::ShowOccupancy()
{
    if(daq_ch)
        SetColor(daq_ch->GetOccupancy());
}

void HyCalModule::ShowEnergy()
{
    ShowEnergy(PRadHyCalModule::GetEnergy());
}

void HyCalModule::ShowEnergy(const double &energy)
{
    SetColor(energy);
}
