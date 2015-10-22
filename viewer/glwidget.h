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

#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include <QtOpenGL>
#include <GL/glu.h>
#include <QDomDocument>
#include <vector>
#include <iostream>
#include <map>
#include <datareader.h>

using namespace std;

enum MainPlain{
    MAIN_PLAIN_XY = 0,
    MAIN_PLAIN_XZ = 1,
    MAIN_PLAIN_YZ = 2
};

struct DisplaySettings{
    double atom_size;
    double max_cutoff;
    double min_cutoff;
    int     atom_quality;
    bool    bHide,bShowDataValues;
    QColor  color;

    DisplaySettings():
        atom_size(1.0),
        atom_quality(5),
        bHide(false),
        bShowDataValues(true),
        max_cutoff(1.0),
        min_cutoff(0.0){}
};

//! [0]
class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget(QWidget *parent = 0);
    ~GLWidget();

    QSize minimumSizeHint() const;
    QSize sizeHint() const;
//! [0]

//! [1]
public slots:
    void setXRotation(double angle);
    void setYRotation(double angle);
    void setXRotation(int angle);
    void setYRotation(int angle);
//    void setZRotation(int angle);

signals:
    void xRotationChanged(int angle);
    void yRotationChanged(int angle);
    void selectedAtoms(int a,int b);
//    void zRotationChanged(int angle);
//! [1]

//! [2]
protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);

//! [2]

//! [3]
private:
    bool bMousePressed,bPickAtom;
    bool bSwapAtomsPressed;

    double xRot;
    double yRot;
    double zRot;
    QPoint lastPos;
    QColor qtPurple;

    QVector3D mass_center;
    QVector3D XY_offset;
    double zoom;
    double ratio;

public:
    int selectedAtomA,selectedAtomB;
    vector<Atom>            *atoms;
    vector<AtomConnection>  *cnts;
    vector<Lead>            *leads;
    AtomData                *data;
    MainPlain mainPlain;
    bool bUseSettingsPerFlag;
    bool bUseOrtho;
    bool bCompileDisplayList;
    bool bSelectedAtomsConnected;
    bool bHideConnections;
    bool bHideAtoms;
    bool bUseLines;
    int  selectedDataColumn,selectedDataSpin,selectDataFlag;


    vector<DisplaySettings> displayPerFlag;
    DisplaySettings          displayAllSettings;
    DisplaySettings          displayConnections;
    map<int,unsigned int> flag2id;
    uint pointTexture;

};
//! [3]

#endif
