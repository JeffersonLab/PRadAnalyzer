#ifndef MARK_SETTING_WIDGET_H
#define MARK_SETTING_WIDGET_H

#include <QWidget>
#include <QColor>
#include <QString>
#include <QStringList>

class QCheckBox;
class QComboBox;
class QLineEdit;
class QSpinBox;
class QPushButton;
class QDoubleSpinBox;

class MarkSettingWidget : public QWidget
{
    Q_OBJECT

public:
    MarkSettingWidget(QWidget *parent = 0);
    MarkSettingWidget(const QStringList &m, QWidget *parent = 0);
    ~MarkSettingWidget() {};

    void SaveSettings();
    void RestoreSettings();

    bool IsChecked();
    bool IsMatchChecked();
    void SetMarkNames(const QStringList &n);
    void AddMarkName(const QString &n);
    void ClearMarkNames();
    void SetChecked(bool c);
    void SetMatchChecked(bool c);
    void SetCurrentMarkIndex(int i);
    void SetWidth(int w);
    void SetColor(const QColor &c);
    void SetColor(const QString &c_text);
    void SetSize(const double &s);

    QString GetCurrentMarkName();
    int GetCurrentMarkIndex();
    int GetWidth();
    QColor GetColor();
    double GetSize();

private slots:
    void colorPicker();

private:
    void initialize();
    void colorLabelChange();

    QComboBox *markCombo;
    QSpinBox *markWidth;
    QDoubleSpinBox *markSize;
    QLineEdit *colorEdit;
    QPushButton *colorButton;
    QCheckBox *checkBox;
    QCheckBox *matchCheckBox;

    QColor chosenColor;
    bool savedCheck;
    bool savedMatch;
    int savedWidth;
    QColor savedColor;
    int savedIndex;
    double savedSize;
};

#endif
