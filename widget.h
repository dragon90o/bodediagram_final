#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>
#include <magnitudeandphase.h>
#include <QtCharts/QChartView>
#include <QPrinter>
#include <QDebug>
#include <QPushButton>
#include <QTextEdit>
#include <QtCharts/QChart>
#include <QStringList>

QT_BEGIN_NAMESPACE
namespace Ui {
class Widget;
}
QT_END_NAMESPACE

class Widget : public QWidget
{
    Q_OBJECT

public:
    Widget(QWidget *parent = nullptr);
    ~Widget();

private slots:

    void onCalculatePushButtonclicked();
    void onExportButtonClicked();

private:
    Ui::Widget *ui;
    MagnitudeAndPhase* mapObject;

    QChart *magnitudeChart = nullptr;
    QChart *phaseChart = nullptr;

};
#endif // WIDGET_H
