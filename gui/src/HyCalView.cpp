//============================================================================//
// A class to contain all the graphic items, including spectrum and modules   //
//                                                                            //
// Chao Peng                                                                  //
// 02/27/2016                                                                 //
//============================================================================//

#include "HyCalView.h"
#include <cmath>
#include <QWheelEvent>
#include <QKeyEvent>

HyCalView::HyCalView(QWidget *parent)
: QGraphicsView(parent)
{
    // enable tracking mode
    setMouseTracking(true);
    viewport()->setMouseTracking(true);
}

void HyCalView::wheelEvent(QWheelEvent *event)
{
    // zoom in / zoom out on mouse wheel scrolling
    double numDegrees = -event->delta() / 8.0;
    double numSteps = numDegrees / 15.0;
    double factor = std::pow(1.125, numSteps);
    scale(factor, factor);
}

void HyCalView::keyPressEvent(QKeyEvent *event)
{
    // press ctrl to enable drag mode
    if(event->key() == Qt::Key_Control) {
        setDragMode(ScrollHandDrag);
    }
}

void HyCalView::keyReleaseEvent(QKeyEvent *event)
{
    // release ctrl to disable drag mode
    if(event->key() == Qt::Key_Control) {
        setDragMode(NoDrag);
    }
}

