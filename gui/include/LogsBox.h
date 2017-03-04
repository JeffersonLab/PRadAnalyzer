#ifndef PRAD_LOG_BOX_H
#define PRAD_LOG_BOX_H

#include <QTextEdit>

#define MAX_NEWLOG_LENGTH 40*256

class QFileSystemWatcher;

class LogsBox : public QTextEdit
{
    Q_OBJECT

public:
    LogsBox(QWidget *parent = NULL);
    virtual ~LogsBox();
    void ShowStdOut();
    void ShowStdErr();
    void TurnOnLog();
    void TurnOffLog();

private slots:
    void handleFileChange(QString path);

private:
    qint64 outpos;
    qint64 errpos;
    QFileSystemWatcher *fileWatcher;
    QString out_path;
    QString err_path;
    FILE *outRedir;
    FILE *errRedir;
    bool logOn;
};

#endif
