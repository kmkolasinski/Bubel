#ifndef DATAREADER_H
#define DATAREADER_H
#include <iostream>
#include <vector>
#include <map>
#include <QVector3D>
#include <QDomDocument>
#include <QFile>
#include <QMessageBox>
#include <QDebug>
#include <string>
#include <sstream>
using namespace std;

struct Atom{
    QVector3D pos;
    int flag,active,no_states,no_bounds;
    Atom():flag(0),active(1),no_states(1),no_bounds(0){
    }
};

struct AtomConnection{
    int atomA,atomB,spinA,spinB;
    double realC,imagC;

    AtomConnection():atomA(0),atomB(0),
                     spinA(1),spinB(1),
                     realC(0),imagC(0){
    }
};

struct Lead{
    bool visible;
    vector<unsigned int> atoms;
    vector<unsigned int> nex_atoms;
    vector<AtomConnection> cnts;
    vector<AtomConnection> inner_cnts;

};

struct AtomsStats{
    QVector3D mass_center;
    double ave_dist;
    QVector<int> flag_list;
    unsigned int no_atoms;
    int max_spin;
    QVector3D min_corner;
    QVector3D max_corner;
    double atom_radius;
    int filter_spin_A,filter_spin_B;
    AtomsStats():mass_center(QVector3D(0,0,0)),ave_dist(0),no_atoms(0),max_spin(0),filter_spin_A(1),filter_spin_B(1)
    {}
};

class DataReader
{
public:
    DataReader();
    ~DataReader();

    void read_data(QString filename);
    void precalculate_data();
private:
    void read_atoms(QDomElement &root);
    void read_connections(QDomElement &root);
    void read_lead(QDomElement &root);
public:
    vector<Atom> atoms;
    vector<AtomConnection> connections;
    vector<Atom> p_atoms;
    vector<AtomConnection> p_connections;
    vector<Lead> leads;
    vector<Lead> p_leads;
    static AtomsStats   atoms_stats;
};

#endif // DATAREADER_H
