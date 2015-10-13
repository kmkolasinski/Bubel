#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    glWidget = new GLWidget;

    ui->horizontalLayoutGL->addWidget(glWidget);

    connect(ui->spinBoxSpinA,SIGNAL(valueChanged(int)),this,SLOT(receiveSpinBoxes(int)));
    connect(ui->spinBoxSpinB,SIGNAL(valueChanged(int)),this,SLOT(receiveSpinBoxes(int)));

    xmlData.read_data("../lattice.xml");
    xmlData.read_data("../lead.xml");
    xmlData.precalculate_data();
    glWidget->atoms = &xmlData.p_atoms;
    glWidget->cnts  = &xmlData.p_connections;
    glWidget->leads = &xmlData.leads;
    update_gui();
}

MainWindow::~MainWindow()
{
    delete glWidget;
    delete ui;
}

QSize MainWindow::sizeHint() const
//! [3] //! [4]
{
    return QSize(1000, 600);
}


void MainWindow::update_gui(){
    QString info;

    for(unsigned i=0; ui->listWidgetLeads->count(); i++){
             QListWidgetItem *item = ui->listWidgetLeads->item(i);
             delete ui->listWidgetLeads->itemWidget(item);
             delete item;
    }

    for(unsigned i=0; i < xmlData.leads.size() ; i++){
            FormLeadParams* lead = new FormLeadParams(i);
            QListWidgetItem *item = new QListWidgetItem();
            item->setSizeHint(QSize(30,40));
            ui->listWidgetLeads->addItem(item);
            ui->listWidgetLeads->setItemWidget(item,lead);
            connect(lead,SIGNAL(emitToggleShowHide(uint,bool)),this,SLOT(toggleLead(uint,bool)));
    }


    AtomsStats& stats = xmlData.atoms_stats;
    info = QString("<b>Number of atoms:</b> ")+QString::number(stats.no_atoms);
    info += QString("<br><b>Flags:</b> ");
    for(unsigned int i = 0 ; i < stats.flag_list.size() ; i++){
        info += QString::number(stats.flag_list[i])+", ";
    }
    info += QString("<br>");
    info += QString("<b>Max spin value:</b> ")+QString::number(stats.max_spin);

    info += QString("<br><b>Ave. dist:</b> ")+QString::number(stats.ave_dist);
    QString str;
    QDebug(&str) << QString("<br><b>Min. pos :</b> ") << stats.min_corner;
    QDebug(&str) << QString("<br><b>Max. pos :</b> ") << stats.max_corner;

    info += str;

    ui->textEditInfo->setText(info);

    ui->spinBoxSpinA->setMinimum(1);
    ui->spinBoxSpinA->setMaximum(stats.max_spin);
    ui->spinBoxSpinB->setMinimum(1);
    ui->spinBoxSpinB->setMaximum(stats.max_spin);

}


void MainWindow::toggleLead(unsigned int id,bool toggle){
    xmlData.leads[id].visible = toggle;
    glWidget->updateGL();

}


void MainWindow::receiveSpinBoxes(int){
    updateWidgets();
}
void MainWindow::updateWidgets(){

    xmlData.atoms_stats.filter_spin_A = ui->spinBoxSpinA->value();
    xmlData.atoms_stats.filter_spin_B = ui->spinBoxSpinB->value();
    xmlData.precalculate_data();
    glWidget->updateGL();
}
