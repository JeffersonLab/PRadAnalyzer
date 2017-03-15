//============================================================================//
// Overload QStyledItemDelegate class to support rich text for QTreeWidget    //
// Based on the example at http://stackoverflow.com/questions/1956542         //
//                                                                            //
// Chao Peng                                                                  //
// 03/02/2016                                                                 //
//============================================================================//

#include "HtmlDelegate.h"

#if QT_VERSION >= 0x050000
#include <QtWidgets>
#else
#include <QtGui>
#endif

void HtmlDelegate::paint(QPainter *painter, const QStyleOptionViewItem &opt, const QModelIndex &index)
const
{
    // self-defined rule, UserRole + 1 means using html delegate or not
    if(!index.data(Qt::UserRole + 1).toBool()) {
        QStyledItemDelegate::paint(painter, opt, index);
        return;
    }

    // copy the option
    QStyleOptionViewItem option = opt;
    initStyleOption(&option, index);

    QStyle *style = option.widget? option.widget->style() : QApplication::style();

    QTextDocument doc;
    doc.setHtml(option.text);

    /// Painting item without text
    option.text = QString();
    style->drawControl(QStyle::CE_ItemViewItem, &option, painter);

    QAbstractTextDocumentLayout::PaintContext ctx;

    // Highlighting text if item is selected
    if (option.state & QStyle::State_Selected)
        ctx.palette.setColor(QPalette::Text, option.palette.color(QPalette::Active, QPalette::HighlightedText));

    QRect textRect = style->subElementRect(QStyle::SE_ItemViewItemText, &option);
    painter->save();
    painter->translate(textRect.topLeft());
    painter->setClipRect(textRect.translated(-textRect.topLeft()));
    doc.documentLayout()->draw(painter, ctx);
    painter->restore();
}

QSize HtmlDelegate::sizeHint(const QStyleOptionViewItem &opt, const QModelIndex &index)
const
{
    QStyleOptionViewItem option = opt;
    initStyleOption(&option, index);

    QTextDocument doc;
    doc.setHtml(option.text);
    doc.setTextWidth(option.rect.width());
    return QSize(doc.idealWidth(), doc.size().height());
}
