#ifndef HTML_DELEGATE_H
#define HTML_DELEGATE_H

#include <QStyledItemDelegate>

#define useHtmlDelegate(i) setData(i, Qt::UserRole + 1, true)

class HtmlDelegate : public QStyledItemDelegate
{
protected:
    void paint ( QPainter * painter, const QStyleOptionViewItem & option, const QModelIndex & index ) const;
    QSize sizeHint ( const QStyleOptionViewItem & option, const QModelIndex & index ) const;
};

#endif
