/****************************************************************************
**
** Copyright (C) 2012 Nokia Corporation and/or its subsidiary(-ies).
** All rights reserved.
** Contact: Nokia Corporation (qt-info@nokia.com)
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
**     the names of its contributors may be used to endorse or promote
**     products derived from this software without specific prior written
**     permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
** $QT_END_LICENSE$
**
****************************************************************************/

#include <QtGui>
#include <QtOpenGL>

#include <math.h>

#include "glwidget.h"


#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE  0x809D
#endif

//! [0]
GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{

    xRot = 0;
    yRot = 0;
    zRot = 0;


    qtPurple = QColor(100,100,100);
    this->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);
    mainPlain = MAIN_PLAIN_XY;
    bUseSettingsPerFlag = false;
    zoom = 1.0;
    bUseOrtho = false;
    bCompileDisplayList = true;

}
//! [0]

//! [1]
GLWidget::~GLWidget()
{

}
//! [1]

//! [2]
QSize GLWidget::minimumSizeHint() const
{
    return QSize(50, 50);
}
//! [2]

//! [3]
QSize GLWidget::sizeHint() const
//! [3] //! [4]
{
    return QSize(600, 600);
}
//! [4]

static void qNormalizeAngle(double &angle)
{
    while (angle < 0)
        angle += 360 ;
    while (angle > 360 )
        angle -= 360 ;
}

//! [5]
void GLWidget::setXRotation(double angle)
{
    qNormalizeAngle(angle);
    if (angle != xRot) {
        xRot = angle;
        emit xRotationChanged(int(angle)%360);
        updateGL();
    }
}
//! [5]

void GLWidget::setYRotation(double angle)
{
    qNormalizeAngle(angle);
    if (angle != yRot) {
        yRot = angle;
        emit yRotationChanged(int(angle)%360);
        updateGL();
    }
}



void GLWidget::setXRotation(int angle)
{
    setXRotation((double)angle);
}
//! [5]

void GLWidget::setYRotation(int angle)
{
    setYRotation((double)angle);
}

//void GLWidget::setZRotation(int angle)
//{
//    qNormalizeAngle(angle);
//    if (angle != zRot) {
//        zRot = angle;
//        emit zRotationChanged(angle);
//        updateGL();
//    }
//}

//! [6]
void GLWidget::initializeGL()
{
    qglClearColor(qtPurple.dark());

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_MULTISAMPLE);
    static GLfloat lightPosition[4] = { 0.5, 5.0, 7.0, 1.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);



}
//! [6]
void renderCylinder(double x1, double y1, double z1, double x2,double y2, double z2, double radius,int subdivisions,GLUquadricObj *quadric)
{

// This is the default direction for the cylinders to face in OpenGL
QVector3D z = QVector3D(0,0,1);
// Get diff between two points you want cylinder along
QVector3D p;
if( qAbs(x1-x2) + qAbs(y1-y2) < 1.0e-4 ){
     p = (QVector3D(x1+0.0001,y1,z1) - QVector3D(x2,y2,z2));
}else p = (QVector3D(x1,y1,z1) - QVector3D(x2,y2,z2));


// Get CROSS product (the axis of rotation)
QVector3D t = QVector3D::crossProduct(z,p);

// Get angle. LENGTH is magnitude of the vector
double angle = 180.0 / M_PI * acos ((QVector3D::dotProduct(z,p) / p.length()));

glPushMatrix();
glTranslated(x2,y2,z2);
glRotated(angle,t.x(),t.y(),t.z());

gluQuadricOrientation(quadric,GLU_OUTSIDE);
gluCylinder(quadric, radius, radius, p.length(), subdivisions, 1);
glPopMatrix();



////handle the degenerate case of z1 == z2 with an approximation
//if(vz == 0)
//    vz = .0001;

//double v = sqrt( vx*vx + vy*vy + vz*vz );
//double ax = 57.2957795*acos( vz/v );
//if ( vz < 0.0 )
//    ax = -ax;
//double rx = -vy*vz;
//double ry = vx*vz;
//glPushMatrix();

////draw the cylinder body
//glTranslatef( x1,y1,z1 );
//glRotatef(ax, rx, ry, 0.0);
//gluQuadricOrientation(quadric,GLU_OUTSIDE);
//gluCylinder(quadric, radius, radius, v, subdivisions, 1);

//glPopMatrix();
}
void renderCylinder_convenient(double x1, double y1, double z1, double x2,double y2, double z2, double radius,int subdivisions)
{
//the same quadric can be re-used for drawing many cylinders
GLUquadricObj *quadric=gluNewQuadric();
gluQuadricNormals(quadric, GLU_SMOOTH);
renderCylinder(x1,y1,z1,x2,y2,z2,radius,subdivisions,quadric);
gluDeleteQuadric(quadric);
}
//! [7]
void GLWidget::paintGL()
{
    static GLuint displayIndex = -1;
    makeCurrent();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double s =zoom;
    if(bUseOrtho)
        glOrtho(-0.5*s*ratio,0.5*s*ratio,-0.5*s,0.5*s, 4.0, 15.0);
    else
        glFrustum(-0.5*s*ratio,0.5*s*ratio,-0.5*s,0.5*s,4.0,15.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    GLfloat d1[4] = { 0.4, 0.5, 0.9, 1.0 };
    GLfloat d2[4] = { 0.9, 0.5, 0.6, 1.0 };
    GLfloat d3[4] = { 0.4, 0.9, 0.2, 1.0 };
    GLfloat d4[4] = { displayAllSettings.color.redF(),
                      displayAllSettings.color.greenF(),
                      displayAllSettings.color.blueF(), 1.0 };
    GLfloat d5[4] = { displayConnections.color.redF(),
                      displayConnections.color.greenF(),
                      displayConnections.color.blueF(), 1.0 };
    GLfloat d6[4] = { 0.8, 0.7, 3.0, 1.0 };
    glMaterialfv(GL_FRONT,GL_DIFFUSE,d4);

    glTranslatef(0.0,0.0,-5.0);
    glRotatef(xRot , 1.0, 0.0, 0.0);
    glRotatef(yRot , 0.0, 1.0, 0.0);
    glRotatef(zRot , 0.0, 0.0, 1.0);


    switch(mainPlain){
        case(MAIN_PLAIN_XY):
        glTranslatef(XY_offset.x(),XY_offset.y(), 0.0);
        break;
        case(MAIN_PLAIN_XZ):
        glTranslatef(XY_offset.y(),0.0,-XY_offset.x());
        break;
        case(MAIN_PLAIN_YZ):
        glTranslatef(0.0,XY_offset.y(),-XY_offset.x());
        break;

    }





    if(bCompileDisplayList){
    glNewList(displayIndex, GL_COMPILE);


    GLUquadricObj *quadric = gluNewQuadric();

    double radius  = displayAllSettings.atom_size * DataReader::atoms_stats.scale * 0.3;
    int no_subdivs = displayAllSettings.atom_quality;

    for(unsigned int i = 0 ; i < atoms->size() ; i++){
        glPushMatrix();
        Atom &atom = (*atoms)[i];
        if(bUseSettingsPerFlag){
            int id = flag2id[atom.flag];
            radius     = displayPerFlag[id].atom_size * DataReader::atoms_stats.scale * 0.3;
            no_subdivs = displayPerFlag[id].atom_quality;
            GLfloat df[4] = { displayPerFlag[id].color.redF(),
                              displayPerFlag[id].color.greenF(),
                              displayPerFlag[id].color.blueF(), 1.0 };
            glMaterialfv(GL_FRONT,GL_DIFFUSE,df);
        }
        glTranslated(atom.pos.x(),atom.pos.y(),atom.pos.z());
        gluSphere( quadric , radius , no_subdivs , no_subdivs );
        glPopMatrix();
    }

    radius     = displayConnections.atom_size * DataReader::atoms_stats.scale * 0.3;
    no_subdivs = displayConnections.atom_quality;

    glMaterialfv(GL_FRONT,GL_DIFFUSE,d5);
    for(unsigned int i = 0 ; i < cnts->size() ; i++){
        AtomConnection &cnt = (*cnts)[i];
        Atom& atomA = (*atoms)[cnt.atomA];
        Atom& atomB = (*atoms)[cnt.atomB];
        renderCylinder(atomA.pos.x(),atomA.pos.y(),atomA.pos.z(),atomB.pos.x(),atomB.pos.y(),atomB.pos.z(),0.2*radius,no_subdivs,quadric);

    }


    glMaterialfv(GL_FRONT,GL_DIFFUSE,d4);

    radius      = displayAllSettings.atom_size * DataReader::atoms_stats.scale * 0.3;
    no_subdivs  = displayAllSettings.atom_quality;

    for(unsigned int l = 0 ; l < leads->size() ; l++){
        Lead& lead = (*leads)[l];
        if(!lead.visible) continue;
        glMaterialfv(GL_FRONT,GL_DIFFUSE,d1);
        for(unsigned int i = 0 ; i < lead.atoms.size() ; i++){
            glPushMatrix();
            Atom &atom = (*atoms)[lead.atoms[i]];
            if(bUseSettingsPerFlag){
                int id = flag2id[atom.flag];
                radius     = displayPerFlag[id].atom_size * DataReader::atoms_stats.scale * 0.3;
                no_subdivs = displayPerFlag[id].atom_quality;
            }

            glTranslated(atom.pos.x(),atom.pos.y(),atom.pos.z());
            gluSphere( quadric , radius*1.1 , 10 , 10 );
            glPopMatrix();
        }
        glMaterialfv(GL_FRONT,GL_DIFFUSE,d2);
        for(unsigned int i = 0 ; i < lead.nex_atoms.size() ; i++){
            glPushMatrix();
            Atom &atom = (*atoms)[lead.nex_atoms[i]];
            if(bUseSettingsPerFlag){
                int id = flag2id[atom.flag];
                radius     = displayPerFlag[id].atom_size * DataReader::atoms_stats.scale * 0.3;
                no_subdivs = displayPerFlag[id].atom_quality;
            }
            glTranslated(atom.pos.x(),atom.pos.y(),atom.pos.z());
            gluSphere( quadric , radius*1.1 , 10 , 10 );
            glPopMatrix();
        }
        glMaterialfv(GL_FRONT,GL_DIFFUSE,d3);
        radius     = displayConnections.atom_size * DataReader::atoms_stats.scale * 0.3;
        no_subdivs = displayConnections.atom_quality;
        for(unsigned int i = 0 ; i < lead.cnts.size() ; i++){
            AtomConnection &cnt = lead.cnts[i];
            Atom& atomA = (*atoms)[cnt.atomA];
            Atom& atomB = (*atoms)[cnt.atomB];
            renderCylinder(atomA.pos.x(),atomA.pos.y(),atomA.pos.z(),atomB.pos.x(),atomB.pos.y(),atomB.pos.z(),0.25*radius,10,quadric);
        }
        glMaterialfv(GL_FRONT,GL_DIFFUSE,d6);
        for(unsigned int i = 0 ; i < lead.inner_cnts.size() ; i++){
            AtomConnection &cnt = lead.inner_cnts[i];
            Atom& atomA = (*atoms)[cnt.atomA];
            Atom& atomB = (*atoms)[cnt.atomB];
            renderCylinder(atomA.pos.x(),atomA.pos.y(),atomA.pos.z(),atomB.pos.x(),atomB.pos.y(),atomB.pos.z(),0.25*radius,10,quadric);
        }
        // this is stupid :/
        lead.shape.indices[0][0] = 0;lead.shape.indices[0][1] = 1;
        lead.shape.indices[1][0] = 1;lead.shape.indices[1][1] = 2;
        lead.shape.indices[2][0] = 2;lead.shape.indices[2][1] = 3;
        lead.shape.indices[3][0] = 3;lead.shape.indices[3][1] = 0;

        lead.shape.indices[4][0] = 4;lead.shape.indices[4][1] = 5;
        lead.shape.indices[5][0] = 5;lead.shape.indices[5][1] = 6;
        lead.shape.indices[6][0] = 6;lead.shape.indices[6][1] = 7;
        lead.shape.indices[7][0] = 7;lead.shape.indices[7][1] = 4;

        lead.shape.indices[8 ][0] = 0;lead.shape.indices[8 ][1] = 4;
        lead.shape.indices[9 ][0] = 1;lead.shape.indices[9 ][1] = 5;
        lead.shape.indices[10][0] = 2;lead.shape.indices[10][1] = 6;
        lead.shape.indices[11][0] = 3;lead.shape.indices[11][1] = 7;

        for(int c = 0 ; c < 12 ; c++){

            int idA = lead.shape.indices[c][0];
            int idB = lead.shape.indices[c][1];

            renderCylinder(lead.shape.data[idA].x(),lead.shape.data[idA].y(),lead.shape.data[idA].z(),
                           lead.shape.data[idB].x(),lead.shape.data[idB].y(),lead.shape.data[idB].z(),
                           0.0015,10,quadric);

        }

    }
    gluDeleteQuadric(quadric);
    bCompileDisplayList  = false;
    glEndList();
    } // end of compile display list;


    glCallList(displayIndex);

}
//! [7]

//! [8]
void GLWidget::resizeGL(int width, int height)
{
    glViewport(0, 0, width, height);
    ratio = double(width)/height;
}
//! [8]

//! [9]
void GLWidget::mousePressEvent(QMouseEvent *event)
{
    lastPos = event->pos();
}
//! [9]

//! [10]
void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    int dx = event->x() - lastPos.x();
    int dy = event->y() - lastPos.y();

    if (event->buttons() & Qt::LeftButton) {
        setXRotation(xRot + 0.2 * dy);
        setYRotation(yRot + 0.2 * dx);
    } else if (event->buttons() & Qt::RightButton) {
        XY_offset += 0.002* QVector3D(dx,-dy,0.0)*zoom;
        updateGL();
    }

    lastPos = event->pos();
}


void GLWidget::wheelEvent(QWheelEvent *event)
{
    float numDegrees = event->delta() / 100.0;

    zoom += numDegrees/20.0;
    zoom = min(max(0.001,zoom),2.0);
    resizeGL(width(),height());
    updateGL();
    event->accept();
}

//! [10]
