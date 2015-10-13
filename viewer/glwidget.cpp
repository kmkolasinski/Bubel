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
#include "qtlogo.h"

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


    qtPurple = QColor::fromCmykF(0.39, 0.39, 0.0, 0.0);
    this->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);
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

static void qNormalizeAngle(int &angle)
{
    while (angle < 0)
        angle += 360 * 16;
    while (angle > 360 * 16)
        angle -= 360 * 16;
}

//! [5]
void GLWidget::setXRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != xRot) {
        xRot = angle;
        emit xRotationChanged(angle);
        updateGL();
    }
}
//! [5]

void GLWidget::setYRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != yRot) {
        yRot = angle;
        emit yRotationChanged(angle);
        updateGL();
    }
}

void GLWidget::setZRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != zRot) {
        zRot = angle;
        emit zRotationChanged(angle);
        updateGL();
    }
}

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
    zoom   = 10.0;

}
//! [6]
void renderCylinder(double x1, double y1, double z1, double x2,double y2, double z2, double radius,int subdivisions,GLUquadricObj *quadric)
{
double vx = x2-x1+0.0000001;
double vy = y2-y1;
double vz = z2-z1;

//handle the degenerate case of z1 == z2 with an approximation
if(vz == 0)
    vz = .0001;

double v = sqrt( vx*vx + vy*vy + vz*vz );
double ax = 57.2957795*acos( vz/v );
if ( vz < 0.0 )
    ax = -ax;
double rx = -vy*vz;
double ry = vx*vz;
glPushMatrix();

//draw the cylinder body
glTranslatef( x1,y1,z1 );
glRotatef(ax, rx, ry, 0.0);
gluQuadricOrientation(quadric,GLU_OUTSIDE);
gluCylinder(quadric, radius, radius, v, subdivisions, 1);

glPopMatrix();
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
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    GLfloat d1[4] = { 0.4, 0.5, 0.9, 1.0 };
    GLfloat d2[4] = { 0.9, 0.5, 0.6, 1.0 };
    GLfloat d3[4] = { 0.4, 0.9, 0.2, 1.0 };
    GLfloat d4[4] = { 0.9, 0.9, 0.7, 1.0 };
    GLfloat d5[4] = { 0.9, 0.9, 1.0, 1.0 };
    GLfloat d6[4] = { 0.8, 0.7, 3.0, 1.0 };
    glMaterialfv(GL_FRONT,GL_DIFFUSE,d4);

    glTranslatef(0.0,0.0,-5.0);
    glRotatef(xRot / 16.0, 1.0, 0.0, 0.0);
    glRotatef(yRot / 16.0, 0.0, 1.0, 0.0);
    glRotatef(zRot / 16.0, 0.0, 0.0, 1.0);
    glTranslatef(XY_offset.x(),XY_offset.y(), 0.0);
    glColor3f( 0.95f, 0.207, 0.031f );
    GLUquadricObj *quadric = gluNewQuadric();
    double radius = DataReader::atoms_stats.atom_radius;
    for(unsigned int i = 0 ; i < atoms->size() ; i++){
        glPushMatrix();
        Atom &atom = (*atoms)[i];
        glTranslated(atom.pos.x(),atom.pos.y(),atom.pos.z());
        gluSphere( quadric , radius , 5 , 5 );
        glPopMatrix();
    }
    glMaterialfv(GL_FRONT,GL_DIFFUSE,d5);
    for(unsigned int i = 0 ; i < cnts->size() ; i++){
        AtomConnection &cnt = (*cnts)[i];
        Atom& atomA = (*atoms)[cnt.atomA];
        Atom& atomB = (*atoms)[cnt.atomB];
        renderCylinder(atomA.pos.x(),atomA.pos.y(),atomA.pos.z(),atomB.pos.x(),atomB.pos.y(),atomB.pos.z(),0.2*radius,4,quadric);

    }


    glMaterialfv(GL_FRONT,GL_DIFFUSE,d4);



    for(unsigned int l = 0 ; l < leads->size() ; l++){
        Lead& lead = (*leads)[l];
        if(!lead.visible) continue;
        glMaterialfv(GL_FRONT,GL_DIFFUSE,d1);
        for(unsigned int i = 0 ; i < lead.atoms.size() ; i++){
            glPushMatrix();
            Atom &atom = (*atoms)[lead.atoms[i]];
            glTranslated(atom.pos.x(),atom.pos.y(),atom.pos.z());
            gluSphere( quadric , radius*1.1 , 10 , 10 );
            glPopMatrix();
        }
        glMaterialfv(GL_FRONT,GL_DIFFUSE,d2);
        for(unsigned int i = 0 ; i < lead.nex_atoms.size() ; i++){
            glPushMatrix();
            Atom &atom = (*atoms)[lead.nex_atoms[i]];
            glTranslated(atom.pos.x(),atom.pos.y(),atom.pos.z());
            gluSphere( quadric , radius*1.1 , 10 , 10 );
            glPopMatrix();
        }
        glMaterialfv(GL_FRONT,GL_DIFFUSE,d3);
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
    }

    delete quadric;
}
//! [7]

//! [8]
void GLWidget::resizeGL(int width, int height)
{
//    int side = qMin(width, height);
    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
//#ifdef QT_OPENGL_ES_1
//    glOrthof(-0.5, +0.5, -0.5, +0.5, 4.0, 15.0);
//#else
    //glOrtho(-0.8, +0.8, -0.8, +0.8, 4.0, 150.0);
//#endif
    gluPerspective(zoom,float(width)/height,0.1,5000);
    glMatrixMode(GL_MODELVIEW);
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
        setXRotation(xRot + 4 * dy);
        setYRotation(yRot + 4 * dx);
    } else if (event->buttons() & Qt::RightButton) {
        XY_offset += 0.005* QVector3D(dx,-dy,0.0)*zoom/20;
        updateGL();
    }

    lastPos = event->pos();
}


void GLWidget::wheelEvent(QWheelEvent *event)
{
    float numDegrees = event->delta() / 100.0;

    zoom += numDegrees;
    zoom = min(max(2.0,zoom),50.0);
    resizeGL(width(),height());
    updateGL();
    event->accept();
}

//! [10]
