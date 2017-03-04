//============================================================================//
// Derived from QTextEdit, redirect stdout and stderr to files                //
// Using file watcher to monitor the files, and display the new input to the  //
// GUI (this file reading way takes resources, but it is currently the only   //
// feasible way I figured out)                                                //
//                                                                            //
// Chao Peng                                                                  //
// 03/13/2016                                                                 //
//============================================================================//

#include <QFileSystemWatcher>
#include <QFile>
#include <QDate>
#include <QTextStream>
#include "LogsBox.h"
#include <cstdio>
#include <unistd.h>

LogsBox::LogsBox(QWidget *parent)
: QTextEdit(parent), fileWatcher(new QFileSystemWatcher(this)), logOn(false)
{
    QDate today = QDate::currentDate();
    QString date = QDate::shortMonthName(today.month()) + tr("_")
                  + QString::number(today.day()) + tr("_")
                  + QString::number(today.year());
    out_path = tr("logs/") + tr("system_") + date + tr(".log");
    err_path = tr("logs/") + tr("error_") + date + tr(".log");

    connect(fileWatcher, SIGNAL(fileChanged(QString)), this, SLOT(handleFileChange(QString)));

    outRedir = freopen(out_path.toStdString().c_str(), "a", stdout);
    setvbuf(stdout,NULL,_IONBF,0);

    errRedir = freopen(err_path.toStdString().c_str(), "a", stderr);
    setvbuf(stderr,NULL,_IONBF,0);

    TurnOnLog();

    QTextEdit::setReadOnly(true);
}


LogsBox::~LogsBox()
{
    stdout = fdopen(STDOUT_FILENO, "w");
    stderr = fdopen(STDERR_FILENO, "w");
    fclose(outRedir);
    fclose(errRedir);
    delete fileWatcher;
}


void LogsBox::TurnOnLog()
{
    fileWatcher->addPath(out_path);
    QFile out_file(out_path);
    outpos = out_file.size();

    fileWatcher->addPath(err_path);
    QFile err_file(err_path);
    errpos = err_file.size();

    logOn = true;
}


void LogsBox::TurnOffLog()
{
    if(!logOn)
        return;

    fileWatcher->removePath(out_path);
    fileWatcher->removePath(err_path);
    logOn = false;
}


void LogsBox::handleFileChange(QString path)
{
    if(path.contains("system"))
        ShowStdOut();

    if(path.contains("error"))
        ShowStdErr();
}


void LogsBox::ShowStdOut()
{
    QFile file(out_path);
    file.open(QFile::ReadOnly | QFile::Text);
    qint64 curpos = file.size();
    if((curpos - outpos) < 0 ||
       ((curpos - outpos) > MAX_NEWLOG_LENGTH))
    {
        QString line = "<font color=\"Black\">Too many log messages, please check system log for more info!</font><br>";
        QTextEdit::moveCursor(QTextCursor::End);
        QTextEdit::insertHtml(line);
        outpos = curpos;
        return;
    }

    file.seek(outpos); 
    outpos = curpos;

    QString header = "<font color=\"Black\">";
    QString end = "</font><br>";

    QTextStream in(&file);

    while(!in.atEnd())
    {
        QString line = in.readLine();
        QTextEdit::moveCursor(QTextCursor::End);
        QTextEdit::insertHtml(header+line+end);
    }
    file.close();
}


void LogsBox::ShowStdErr()
{
    QFile file(err_path);
    file.open(QFile::ReadOnly | QFile::Text);
    qint64 curpos = file.size();

    if((curpos - errpos) < 0 ||
       ((curpos - errpos) > MAX_NEWLOG_LENGTH))
    {
        QString line = "<font color=\"Red\">Too many error messages, please check error log for more info!</font><br>";
        QTextEdit::moveCursor(QTextCursor::End);
        QTextEdit::insertHtml(line);
        errpos = curpos;
        return;
    }

    file.seek(errpos);
    errpos = curpos;

    QString header = "<font color=\"Red\">";
    QString end = "</font><br>";

    QTextStream in(&file);

    while(!in.atEnd())
    {
        QString line = in.readLine();
        QTextEdit::moveCursor(QTextCursor::End);
        QTextEdit::insertHtml(header+line+end);
    }
    file.close();
}
