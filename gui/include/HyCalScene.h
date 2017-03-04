#ifndef QT_HYCAL_SCENE_H
#define QT_HYCAL_SCENE_H

#include <QGraphicsScene>
#include <vector>
#include "PRadHyCalDetector.h"

class PRadEventViewer;
class HyCalModule;

class HyCalScene : public QGraphicsScene, public PRadHyCalDetector
{
    Q_OBJECT

public:
    struct TextBox
    {
        QString name;
        QString text;
        QColor textColor;
        QRectF bound;
        QColor bgColor;

        TextBox() {};
        TextBox(const QString &n, const QString &t, const QColor &tc, const QRectF &b, const QColor &c)
        : name(n), text(t), textColor(tc), bound(b), bgColor(c) {};
        TextBox(const QString &n, const QColor &tc, const QRectF &b, const QColor &c)
        : name(n), text(""), textColor(tc), bound(b), bgColor(c) {};
    };

    struct MarkAttributes
    {
        int shape_index;
        int width;
        QColor color;
        double size;

        MarkAttributes()
        : shape_index(0), width(1), color(Qt::black), size(5.0)
        {};
        MarkAttributes(const int &i, const int &w, const QColor &c, const double &s)
        : shape_index(i), width(w), color(c), size(s)
        {};
    };
    static QStringList GetShapeList();

    struct HitsMark
    {
        QString name;
        QString text;
        QColor textColor;
        QRectF textBox;
        QPointF hitPos;
        MarkAttributes attr;

        HitsMark() {};
        HitsMark(const QString &n, const QString &t, const QPointF &p, const MarkAttributes& m)
        : name(n), text(t), textColor(Qt::black), textBox(QRectF(p.x()-50., p.y()-20., 100., 40.)),
          hitPos(p), attr(m)
        {};
    };

public:
    HyCalScene(PRadEventViewer *p, QObject *parent = 0)
    : QGraphicsScene(parent), console(p),
      pModule(nullptr), sModule(nullptr), rModule(nullptr), showScalers(false)
    {};
    HyCalScene(PRadEventViewer *p, const QRectF & sceneRect, QObject *parent = 0)
    : QGraphicsScene(sceneRect, parent), console(p),
      pModule(nullptr), sModule(nullptr), rModule(nullptr), showScalers(false)
    {};
    HyCalScene(PRadEventViewer*p, qreal x, qreal y, qreal width, qreal height, QObject *parent = 0)
    : QGraphicsScene(x, y, width, height, parent), console(p),
      pModule(nullptr), sModule(nullptr), rModule(nullptr), showScalers(false)
    {};

    void AddTDCBox(const QString &name,
                   const QColor &textColor,
                   const QRectF &textBox,
                   const QColor &bgColor);
    void AddScalerBox(const QString &name,
                      const QColor &textColor,
                      const QRectF &textBox,
                      const QColor &bgColor);
    void AddHitsMark(const QString &name,
                     const QPointF &position,
                     const MarkAttributes &m = MarkAttributes(),
                     const QString &text = "");
    void ClearHitsMarks();
    void UpdateScalerBox(const QString &text, const int &group = 0);
    void UpdateScalerBox(const QStringList &texts);
    void ShowScalers(const bool &s = true) {showScalers = s;};
    void ShowCluster(const ModuleCluster &cluster);
    void ShowCluster(int index);
    template<typename... Args>
    void ModuleAction(void (HyCalModule::*act)(Args...), Args&&... args)
    {
        for(auto &module : module_list)
            (((HyCalModule*)module)->*act)(std::forward<Args>(args)...);
    }
    // overloaded
    void ReadModuleList(const std::string &path);

protected:
    void drawForeground(QPainter *painter, const QRectF &rect);
    void mousePressEvent(QGraphicsSceneMouseEvent *event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);

private:
    void drawScalerBoxes(QPainter *painter);
    void drawTDCBoxes(QPainter *painter);
    void drawHitsMarks(QPainter *painter);
    void drawHitsMark(QPainter *painter, const QPointF &pos, const MarkAttributes &attr);

private:
    PRadEventViewer *console;
    HyCalModule *pModule, *sModule, *rModule;
    bool showScalers;
    QList<TextBox> tdcBoxList;
    QVector<TextBox> scalarBoxList;
    QVector<HitsMark> hitsMarkList;
};

#endif
