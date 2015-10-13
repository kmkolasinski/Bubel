#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <glwidget.h>
#include <formleadparams.h>
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    QSize sizeHint() const;
    ~MainWindow();
public slots:
    void receiveSpinBoxes(int);
    void receiveDoubleSpinBoxes(double);
    void receiveCheckBoxes(int);

    void updateWidgets();
    void toggleLead(unsigned int id,bool toggle);
    void open();
    void chooseColor();
    void chooseConnectionColor();
    void togglePerFlagDisplaySettings(bool);
    void updatePerFlagSettings(int);

private:
    bool bSkipSignals;
    void update_gui();

    Ui::MainWindow *ui;
    GLWidget* glWidget;
    DataReader xmlData;

};

#endif // MAINWINDOW_H
