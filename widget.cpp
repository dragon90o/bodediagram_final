#include "widget.h"
#include "ui_widget.h"
#include "magnitudeandphase.h"

Widget::Widget(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::Widget)
    , mapObject(ui->transferFunctionlineEdit, ui->sigmaLineEdit, ui->frequencyMaxlineEdit, ui->frequencyMinlineEdit, this)
{
    ui->setupUi(this);
    connect(ui->calculatePushButton, &QPushButton::clicked, this, &Widget::onCalculatePushButtonclicked);
}

Widget::~Widget()
{
    delete ui;
}

void Widget::onCalculatePushButtonclicked()
{
    mapObject.frequencies();
}
