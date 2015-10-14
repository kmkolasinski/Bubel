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
    connect(ui->pushButtonOpen,SIGNAL(released()),this,SLOT(open()));

    QButtonGroup* group = new QButtonGroup(this);
    group->addButton(ui->radioButtonMainXY);
    group->addButton(ui->radioButtonMainXZ);
    group->addButton(ui->radioButtonMainYZ);
    connect(group,SIGNAL(buttonClicked(int)),this,SLOT(receiveCheckBoxes(int)));
    connect(ui->checkBoxOrthoProj,SIGNAL(stateChanged(int)),this,SLOT(receiveCheckBoxes(int)));

    connect(ui->doubleSpinBoxAtomSize,SIGNAL(valueChanged(double)),this,SLOT(receiveDoubleSpinBoxes(double)));
    connect(ui->spinBoxAtomQ,SIGNAL(valueChanged(int)),this,SLOT(receiveSpinBoxes(int)));
    connect(ui->pushButtonAtomColor,SIGNAL(released()),this,SLOT(chooseColor()));

    connect(ui->doubleSpinBoxConnectionSize,SIGNAL(valueChanged(double)),this,SLOT(receiveDoubleSpinBoxes(double)));
    connect(ui->spinBoxConnectionQ,SIGNAL(valueChanged(int)),this,SLOT(receiveSpinBoxes(int)));
    connect(ui->pushButtonConnectionColour,SIGNAL(released()),this,SLOT(chooseConnectionColor()));


    connect(ui->checkBoxPerFlagSettings,SIGNAL(toggled(bool)),this,SLOT(togglePerFlagDisplaySettings(bool)));

    connect(ui->comboBoxFlags,SIGNAL(currentIndexChanged(int)),this,SLOT(updatePerFlagSettings(int)));


    connect(glWidget,SIGNAL(xRotationChanged(int)),ui->spinBoxRotX,SLOT(setValue(int)));
    connect(glWidget,SIGNAL(yRotationChanged(int)),ui->spinBoxRotY,SLOT(setValue(int)));


    connect(ui->spinBoxRotX,SIGNAL(valueChanged(int)),glWidget,SLOT(setXRotation(int)));
    connect(ui->spinBoxRotY,SIGNAL(valueChanged(int)),glWidget,SLOT(setYRotation(int)));

    xmlData.read_data("../lattice.xml");
    xmlData.read_data("../lead.xml");
    xmlData.precalculate_data();
    glWidget->atoms = &xmlData.p_atoms;
    glWidget->cnts  = &xmlData.p_connections;
    glWidget->leads = &xmlData.p_leads;
    bSkipSignals = false;
    lastDir = "";
    update_gui();

    glWidget->setAcceptDrops(true);
    setAcceptDrops(true);

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

    bSkipSignals = true;
    ui->doubleSpinBoxAtomSize->setValue(xmlData.atoms_stats.atom_radius);
    ui->doubleSpinBoxConnectionSize->setValue(xmlData.atoms_stats.atom_radius);

    ui->doubleSpinBoxAtomSize->setSingleStep(xmlData.atoms_stats.atom_radius/10.0);
    ui->doubleSpinBoxConnectionSize->setSingleStep(xmlData.atoms_stats.atom_radius/10.0);

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
    glWidget->displayPerFlag.clear();
    glWidget->flag2id.clear();
    ui->comboBoxFlags->clear();
    for(unsigned int i = 0 ; i < stats.flag_list.size() ; i++){
        info += QString::number(stats.flag_list[i])+", ";
        glWidget->flag2id[stats.flag_list[i]] = i;
        DisplaySettings ds;
        ds.atom_quality  = ui->spinBoxAtomQ->value();
        ds.atom_size     = ui->doubleSpinBoxAtomSize->value();
        QPalette palette = ui->pushButtonAtomColor->palette();
        QColor color = palette.color(QPalette::Button);
        ds.color = color;

        glWidget->displayPerFlag.push_back(ds);
        ui->comboBoxFlags->addItem("Flag="+QString::number(stats.flag_list[i]));
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

    bSkipSignals = false;
    updateWidgets();
}

void MainWindow::toggleLead(unsigned int id,bool toggle){
    xmlData.leads[id].visible = toggle;
    updateWidgets();
}

void MainWindow::receiveSpinBoxes(int){
    updateWidgets();
}

void MainWindow::receiveDoubleSpinBoxes(double){
    updateWidgets();
}

void MainWindow::receiveCheckBoxes(int toggle){
    toggle = 0;
    if(ui->radioButtonMainXY->isChecked()) glWidget->mainPlain = MAIN_PLAIN_XY;
    if(ui->radioButtonMainXZ->isChecked()) glWidget->mainPlain = MAIN_PLAIN_XZ;
    if(ui->radioButtonMainYZ->isChecked()) glWidget->mainPlain = MAIN_PLAIN_YZ;

    glWidget->bUseOrtho = ui->checkBoxOrthoProj->isChecked();
    glWidget->updateGL();
}

void MainWindow::updateWidgets(){
    if(bSkipSignals) return;
    glWidget->bCompileDisplayList = true;
    xmlData.atoms_stats.filter_spin_A = ui->spinBoxSpinA->value();
    xmlData.atoms_stats.filter_spin_B = ui->spinBoxSpinB->value();
    xmlData.precalculate_data();

    glWidget->displayConnections.atom_quality = ui->spinBoxConnectionQ->value();
    glWidget->displayConnections.atom_size    = ui->doubleSpinBoxConnectionSize->value();
    QPalette palette = ui->pushButtonConnectionColour->palette();
    QColor color = palette.color(QPalette::Button);
    glWidget->displayConnections.color = color;


    if(!ui->checkBoxPerFlagSettings->isChecked()){
        glWidget->displayAllSettings.atom_quality = ui->spinBoxAtomQ->value();
        glWidget->displayAllSettings.atom_size    = ui->doubleSpinBoxAtomSize->value();
        QPalette palette = ui->pushButtonAtomColor->palette();
        QColor color = palette.color(QPalette::Button);
        glWidget->displayAllSettings.color = color;
    }else{
        int index = ui->comboBoxFlags->currentIndex();
        glWidget->displayPerFlag[index].atom_quality = ui->spinBoxAtomQ->value();
        glWidget->displayPerFlag[index].atom_size    = ui->doubleSpinBoxAtomSize->value();
        QPalette palette = ui->pushButtonAtomColor->palette();
        QColor color = palette.color(QPalette::Button);
        glWidget->displayPerFlag[index].color = color;
    }
    glWidget->updateGL();
}

void MainWindow::open(){
    QString fn = QFileDialog::getOpenFileName(this, tr("Open Quantulaba XML File"),
                                                  lastDir, tr("xml (*.xml);;All Files (*)"));
    if (!fn.isEmpty()){
        lastDir = fn;
        xmlData.read_data(fn);
        xmlData.precalculate_data();
        update_gui();
    }
}

void MainWindow::dragEnterEvent(QDragEnterEvent *event)
{

    if (event->mimeData()->hasFormat("text/plain"))
        event->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent *event)
{
    QList<QUrl> urls = event->mimeData()->urls();
    if (urls.isEmpty())
        return;

    for(int i = 0 ; i < urls.size() ; i++){
        QString fileName = urls[i].toLocalFile();
        QFileInfo fileInfo(fileName);
        if(fileInfo.suffix() == "xml"){
            qDebug() << "Opening:" << fileName ;
            lastDir = fileName;
            xmlData.read_data(fileName);
        }
    }
    xmlData.precalculate_data();
    update_gui();

}


void MainWindow::chooseColor(){


    QPalette palette = ui->pushButtonAtomColor->palette();
    QColor color = palette.color(QPalette::Button);

    QColorDialog colorDialog;
    colorDialog.setCurrentColor(color);
    colorDialog.setParent(this);
    if(colorDialog.exec()){
        color = colorDialog.selectedColor();
        palette.setColor(QPalette::Button, color);
        ui->pushButtonAtomColor->setAutoFillBackground(true);
        ui->pushButtonAtomColor->setPalette(palette);
        updateWidgets();
    }
}
void MainWindow::chooseConnectionColor(){
    QPalette palette = ui->pushButtonConnectionColour->palette();
    QColor color = palette.color(QPalette::Button);

    QColorDialog colorDialog;
    colorDialog.setCurrentColor(color);
    colorDialog.setParent(this);
    if(colorDialog.exec()){
        color = colorDialog.selectedColor();
        palette.setColor(QPalette::Button, color);
        ui->pushButtonConnectionColour->setAutoFillBackground(true);
        ui->pushButtonConnectionColour->setPalette(palette);
        updateWidgets();
    }
}

void MainWindow::togglePerFlagDisplaySettings(bool toggle){

    ui->comboBoxFlags->setEnabled(toggle);
    glWidget->bUseSettingsPerFlag = toggle;
    updateWidgets();

}

void MainWindow::updatePerFlagSettings(int row){
    if(bSkipSignals){
        return;
    }
    bSkipSignals = true;
    ui->doubleSpinBoxAtomSize->setValue(glWidget->displayPerFlag[row].atom_size);
    ui->spinBoxAtomQ->setValue(glWidget->displayPerFlag[row].atom_quality);
    QPalette palette;
    palette.setColor(QPalette::Button, glWidget->displayPerFlag[row].color);
    ui->pushButtonAtomColor->setPalette(palette);

    bSkipSignals = false;
    updateWidgets();
}
