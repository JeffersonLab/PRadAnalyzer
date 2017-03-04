//============================================================================//
// A class derived from QGraphicsScene                                        //
// This class is used for self-defined selection behavior                     //
// It also displays the reconstructed clusters                                //
//                                                                            //
// Chao Peng, Weizhi Xiong                                                    //
// 02/27/2016                                                                 //
//============================================================================//

#include "HyCalScene.h"
#include "PRadEventViewer.h"
#include "HyCalModule.h"

#if QT_VERSION >= 0x050000
#include <QtWidgets>
#else
#include <QtGui>
#endif

// static function to return the supported mark shapes
QStringList HyCalScene::GetShapeList()
{
    QStringList list;
    list << "Circle" << "Diamond" << "Square" << "Triangle" << "Cross"
         << "X-Cross";
    return list;
}

void HyCalScene::drawForeground(QPainter *painter, const QRectF &rect)
{
    QGraphicsScene::drawForeground(painter, rect);

    painter->save();

    if(showScalers)
        drawScalerBoxes(painter);

    if(console->GetAnnoType() == ShowTDC)
        drawTDCBoxes(painter);

    if(console->GetViewMode() == EnergyView)
        drawHitsMarks(painter);

    painter->restore();
}


// show scaler boxes
void HyCalScene::drawScalerBoxes(QPainter *painter)
{
    painter->setFont(QFont("times", 16, QFont::Bold));
    for(auto it = scalarBoxList.begin(); it != scalarBoxList.end(); ++it)
    {
        QPen pen(it->textColor);
        pen.setWidth(2);
        pen.setCosmetic(true);
        painter->setPen(pen);
        QPainterPath path;
        path.addRect(it->bound);
        painter->fillPath(path, it->bgColor);
        painter->drawPath(path);
        QRectF name_box = it->bound.translated(0, -it->bound.height());
        painter->drawText(name_box,
                          it->name,
                          QTextOption(Qt::AlignBottom | Qt::AlignHCenter));
        painter->drawText(it->bound,
                          it->text,
                          QTextOption(Qt::AlignCenter | Qt::AlignHCenter));
    }
}

// show the tdc groups
void HyCalScene::drawTDCBoxes(QPainter *painter)
{
    painter->setFont(QFont("times", 24, QFont::Bold));
    for(auto it = tdcBoxList.begin(); it != tdcBoxList.end(); ++it)
    {
        QPen pen(it->textColor);
        pen.setWidth(2);
        pen.setCosmetic(true);
        painter->setPen(pen);
        QPainterPath path;
        path.addRect(it->bound);
        painter->fillPath(path, it->bgColor);
        painter->drawPath(path);
        painter->drawText(it->bound,
                          it->name,
                          QTextOption(Qt::AlignCenter | Qt::AlignHCenter));

    }
}

// draw hits marks
void HyCalScene::drawHitsMarks(QPainter *painter)
{
    for(auto &mark : hitsMarkList)
    {
        drawHitsMark(painter, mark.hitPos, mark.attr);

        // draw text
        if(mark.text.isEmpty())
            continue;

        QPen pen(mark.textColor);
        painter->setPen(pen);
        painter->setFont(QFont("times", 9));
        painter->drawText(mark.textBox, mark.text, QTextOption(Qt::AlignTop | Qt::AlignHCenter));
    }
}

// helper, draw a single mark
void HyCalScene::drawHitsMark(QPainter *painter, const QPointF& pos, const HyCalScene::MarkAttributes &attr)
{
    // draw a circle mark
    QPen pen(attr.color);
    pen.setWidth(attr.width);
    painter->setPen(pen);

    switch(attr.shape_index)
    {
    default: // not supported
        break;
    case 0: // circle
        painter->drawEllipse(pos, attr.size, attr.size);
        break;
    case 1: // diamond
      {
        QPointF sqr[4] = {
            QPointF(pos.x() + attr.size, pos.y()),
            QPointF(pos.x(), pos.y() + attr.size),
            QPointF(pos.x() - attr.size, pos.y()),
            QPointF(pos.x(), pos.y() - attr.size)
        };
        painter->drawConvexPolygon(sqr, 4);
      }
        break;
    case 2: // square
      {
        QPointF sqr[4] = {
            QPointF(pos.x() + attr.size, pos.y() - attr.size),
            QPointF(pos.x() - attr.size, pos.y() - attr.size),
            QPointF(pos.x() - attr.size, pos.y() + attr.size),
            QPointF(pos.x() + attr.size, pos.y() + attr.size)
        };
        painter->drawConvexPolygon(sqr, 4);
      }
        break;
    case 3: // triangle
      {
        QPointF tri[3] = {
            QPointF(pos.x(), pos.y() - attr.size),
            QPointF(pos.x() - attr.size, pos.y() + attr.size),
            QPointF(pos.x() + attr.size, pos.y() + attr.size)
        };
        painter->drawConvexPolygon(tri, 3);
      }
        break;
    case 4: // cross
      {
        QPainterPath cross;
        cross.moveTo(pos.x(), pos.y() - attr.size);
        cross.lineTo(pos.x(), pos.y() + attr.size);
        cross.moveTo(pos.x() + attr.size, pos.y());
        cross.lineTo(pos.x() - attr.size, pos.y());
        painter->drawPath(cross);
      }
        break;
    case 5: // X-cross
      {
        QPainterPath cross;
        cross.moveTo(pos.x() - attr.size, pos.y() - attr.size);
        cross.lineTo(pos.x() + attr.size, pos.y() + attr.size);
        cross.moveTo(pos.x() + attr.size, pos.y() - attr.size);
        cross.lineTo(pos.x() - attr.size, pos.y() + attr.size);
        painter->drawPath(cross);
      }
        break;
    }
}

// pverloaded read module list
void HyCalScene::ReadModuleList(const std::string &path)
{
    if(path.empty())
        return;

    ConfigParser c_parser;
    if(!c_parser.ReadFile(path)) {
        std::cerr << "PRad HyCal Detector Error: Failed to read module list file "
                  << "\"" << path << "\"."
                  << std::endl;
        return;
    }

    // clear all modules
    ClearModuleList();

    std::string name;
    std::string type, sector;
    PRadHyCalModule::Geometry geo;

    // some info that is not read from list
    while (c_parser.ParseLine())
    {
        if(!c_parser.CheckElements(8))
            continue;

        c_parser >> name >> type
                 >> geo.size_x >> geo.size_y >> geo.size_z
                 >> geo.x >> geo.y >> geo.z;

        geo.type = PRadHyCalModule::get_module_type(type.c_str());

        HyCalModule *module = new HyCalModule(console, name, geo);

        // failed to add module to detector
        if(!AddModule(module)) {
            delete module;
            continue;
        }
        addItem(module);
    }

    // sort the module by id
    SortModuleList();
}


void HyCalScene::AddTDCBox(const QString &text,
                           const QColor &textColor,
                           const QRectF &textBox,
                           const QColor &bkgColor)
{
    tdcBoxList.append(TextBox(text, textColor, textBox, bkgColor));
}

void HyCalScene::AddScalerBox(const QString &text,
                              const QColor &textColor,
                              const QRectF &textBox,
                              const QColor &bkgColor)
{
    scalarBoxList.push_back(TextBox(text, textColor, textBox, bkgColor));
}

void HyCalScene::AddHitsMark(const QString &name,
                             const QPointF &position,
                             const HyCalScene::MarkAttributes &m,
                             const QString &text)
{
    hitsMarkList.push_back(HitsMark(name, text, position, m));
}

void HyCalScene::UpdateScalerBox(const QString &text, const int &group)
{
    if(group < 0 || group >= scalarBoxList.size())
        return;

    scalarBoxList[group].text = text;
}

void HyCalScene::UpdateScalerBox(const QStringList &texts)
{
    for(int i = 0; i < texts.size(); ++i)
    {
        UpdateScalerBox(texts.at(i), i);
    }
}

// overload to change the default selection behavior
void HyCalScene::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
    if(event->modifiers().testFlag(Qt::ControlModifier))
        return;
//    QGraphicsScene::mousePressEvent(event);
    if(event->button() == Qt::LeftButton)
        pModule = dynamic_cast<HyCalModule*>(itemAt(event->scenePos(), QTransform()));
    else if(event->button() == Qt::RightButton)
        rModule = dynamic_cast<HyCalModule*>(itemAt(event->scenePos(), QTransform()));
}

// overload to change the default selection behavior
void HyCalScene::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
    if(event->button() == Qt::LeftButton) {
        if(pModule &&
           !pModule->isSelected() &&
           pModule == dynamic_cast<HyCalModule*>(itemAt(event->scenePos(), QTransform()))
          ) {
            pModule->setSelected(true);
            if(sModule) {
                sModule->setSelected(false);
            }
            sModule = pModule;
        }
    } else if(event->button() == Qt::RightButton) {
        if(rModule &&
           rModule->isSelected() &&
           rModule == dynamic_cast<HyCalModule*>(itemAt(event->scenePos(), QTransform()))
          ) {
            rModule->setSelected(false);
            sModule = nullptr;
        }
    }
}

void HyCalScene::ClearHitsMarks()
{
    hitsMarkList.clear();
}

void HyCalScene::ShowCluster(int index)
{
    if((size_t)index >= module_clusters.size())
        ModuleAction(&HyCalModule::ShowEnergy);

    ShowCluster(module_clusters.at(index));
}

void HyCalScene::ShowCluster(const ModuleCluster &cluster)
{
    // erase all the modules
    for(auto module : module_list)
    {
        ((HyCalModule*)module)->ShowEnergy(0);
    }

    // show only the cluster info
    for(auto hit : cluster.hits)
    {
        HyCalModule *module = (HyCalModule*)PRadHyCalDetector::GetModule(hit.id);
        if(module)
            module->ShowEnergy(hit.energy);
    }
}

