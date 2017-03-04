#ifndef QT_HYCAL_VIEW_H
#define QT_HYCAL_VIEW_H

#include <QGraphicsView>

class HyCalView : public QGraphicsView
{
    Q_OBJECT

public:
    HyCalView(QWidget *parent = 0);

protected:
    void wheelEvent(QWheelEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void keyReleaseEvent(QKeyEvent *event);
};

#endif
