#include "datareader.h"
AtomsStats  DataReader::atoms_stats;

DataReader::DataReader()
{
    atoms.clear();
    p_atoms.clear();
    connections.clear();
    p_connections.clear();
}

DataReader::~DataReader()
{
    atoms.clear();
    p_atoms.clear();
    connections.clear();
    p_connections.clear();
}

void DataReader::read_data(QString filename){


    QFile* xmlFile = new QFile(filename);
            if (!xmlFile->open(QIODevice::ReadOnly | QIODevice::Text)) {
                    QMessageBox::critical(NULL,"Load XML File Problem",
                         "Couldn't open xml to load settings for download",
                    QMessageBox::Ok);
                    return;
            }
    QDomDocument xmlBOM;
    xmlBOM.setContent(xmlFile);

    QDomElement root=xmlBOM.documentElement();
    if(root.tagName() == "lattice") {
        qDebug() << "#Reading lattice:";
        QDomElement Component=root.firstChild().toElement();
        while(!Component.isNull())
        {
            // Check if the child tag name is COMPONENT
            if (Component.tagName()=="atoms")
            {
                read_atoms(Component);
            }
            // Check if the child tag name is COMPONENT
            if (Component.tagName()=="connections")
            {
                read_connections(Component);
            }

            Component = Component.nextSibling().toElement();
        }
    }else if(root.tagName() == "atoms") read_atoms(root);
    else if(root.tagName() == "lead"){
        leads.clear();
        read_lead(root);
    }
    xmlFile->close();

}


void DataReader::read_atoms(QDomElement& root){
    atoms.clear();
    atoms_stats.flag_list.clear();
    atoms_stats.max_spin = 0;




    qDebug() << "Reading atoms:";
    qDebug() << root.tagName();
    QDomElement Component=root.firstChild().toElement();

    // Loop while there is a child
    while(!Component.isNull())
    {
        // Check if the child tag name is COMPONENT
        if (Component.tagName()=="d")
        {
            double position[3];
            int flag,active,no_states,no_bounds;
            // Read each child of the component node
            QString line;
            line=Component.firstChild().toText().data();


            std::stringstream stream(line.toStdString());
            stream >> position[0];
            stream >> position[1];
            stream >> position[2];
            stream >> flag;
            stream >> active;
            stream >> no_states;
            stream >> no_bounds;
            Atom atom;
            atom.pos       = QVector3D(position[0],position[1],position[2]);
            atom.flag      = flag;
            atom.active    = active;
            atom.no_states = no_states;
            atom.no_bounds = no_bounds;

            atoms.push_back(atom);
            atoms_stats.mass_center += atom.pos;

            if(!atoms_stats.flag_list.contains(flag)) atoms_stats.flag_list.push_back(flag);
            if(atoms_stats.max_spin < no_states) atoms_stats.max_spin = no_states;
            atoms_stats.min_corner = atom.pos;
            atoms_stats.max_corner = atom.pos;



        }
        // Next component
        Component = Component.nextSibling().toElement();
    }

    for(unsigned int i = 0 ; i < atoms.size() ; i++){
        // bbox
        if(atoms_stats.min_corner.x() > atoms[i].pos.x()) atoms_stats.min_corner.setX(atoms[i].pos.x());
        if(atoms_stats.min_corner.y() > atoms[i].pos.y()) atoms_stats.min_corner.setY(atoms[i].pos.y());
        if(atoms_stats.min_corner.z() > atoms[i].pos.z()) atoms_stats.min_corner.setZ(atoms[i].pos.z());

        if(atoms_stats.max_corner.x() < atoms[i].pos.x()) atoms_stats.max_corner.setX(atoms[i].pos.x());
        if(atoms_stats.max_corner.y() < atoms[i].pos.y()) atoms_stats.max_corner.setY(atoms[i].pos.y());
        if(atoms_stats.max_corner.z() < atoms[i].pos.z()) atoms_stats.max_corner.setZ(atoms[i].pos.z());
    }
    atoms_stats.mass_center  /= atoms.size();
    atoms_stats.no_atoms      = atoms.size();
}

void DataReader::read_connections(QDomElement& root){
    connections.clear();

    qDebug() << "Reading connections:";
    qDebug() << root.tagName();
    QDomElement Component=root.firstChild().toElement();
    QVector3D dist;
    int iter = 0;

    atoms_stats.ave_dist = 0;

    // Loop while there is a child
    while(!Component.isNull())
    {
        // Check if the child tag name is COMPONENT
        if (Component.tagName()=="d")
        {
            // Get the first child of the component
            //QDomElement Child=Component.firstChild().toElement();

            AtomConnection connection;

            QString line = Component.firstChild().toText().data();
            std::stringstream stream(line.toStdString());
            stream >> connection.atomA;
            stream >> connection.atomB;
            stream >> connection.spinA;
            stream >> connection.spinB;
            stream >> connection.realC;
            stream >> connection.imagC;

            connection.atomA--;
            connection.atomB--;

//            // Read each child of the component node
//            while (!Child.isNull())
//            {
//                // Read Name and value
//                if (Child.tagName()=="A")   connection.atomA=Child.firstChild().toText().data().toInt()-1;
//                if (Child.tagName()=="B")   connection.atomB=Child.firstChild().toText().data().toInt()-1;
//                if (Child.tagName()=="sA")  connection.spinA=Child.firstChild().toText().data().toInt();
//                if (Child.tagName()=="sB")  connection.spinB=Child.firstChild().toText().data().toInt();
//                if (Child.tagName()=="vr")  connection.realC=Child.firstChild().toText().data().toDouble();
//                if (Child.tagName()=="vi")  connection.imagC=Child.firstChild().toText().data().toDouble();
//                // Next child
//                Child = Child.nextSibling().toElement();
//            }

            connections.push_back(connection);
            if(connection.atomA != connection.atomB){
                dist = atoms[connection.atomA].pos - atoms[connection.atomB].pos;
                atoms_stats.ave_dist += dist.length();
                iter++;
            }
        }
        // Next component
        Component = Component.nextSibling().toElement();
    }
    atoms_stats.ave_dist /= iter;
}


void DataReader::read_lead(QDomElement &root){

    qDebug() << "Reading lead:";
    qDebug() << root.tagName();
    QDomElement Component=root.firstChild().toElement();

    Lead lead;
    bool shapeTypeInitialized = false;


    // Loop while there is a child
    while(!Component.isNull())
    {

        // Check if the child tag name is COMPONENT
        if (Component.tagName()=="shape_type")
        {


            qDebug() << "Lead area:" << Component.firstChild().toText().data();
            if(Component.firstChild().toText().data()      == "SHAPE_RECTANGLE_XY")     lead.shape.type = SHAPE_RECTANGLE_XY;
            else if(Component.firstChild().toText().data() == "SHAPE_CONVEX_QUAD_XY")   lead.shape.type = SHAPE_CONVEX_QUAD_XY;
            else if(Component.firstChild().toText().data() == "SHAPE_RANGE_3D")         lead.shape.type = SHAPE_RANGE_3D;
            else{
                lead.shape.type = SHAPE_NONE;
                qDebug() << "Warning: Shape:" << Component.firstChild().toText().data() << " is not defined." ;
            }

            shapeTypeInitialized = true;

        }


        // Check if the child tag name is COMPONENT
        if (Component.tagName()=="shape_data")
        {

            qDebug() << "Reading shape data";
            if(!shapeTypeInitialized){
                qDebug() << "Warning::The shape_type keyword must be before shape_data in XML file.";
                return;
            }

            QString line;
            line=Component.firstChild().toText().data();
            std::stringstream stream(line.toStdString());
            double corners[8][3] = {0};

            switch(lead.shape.type){
                case(SHAPE_RECTANGLE_XY):
                    stream >> corners[0][0];
                    stream >> corners[0][1];
                    stream >> corners[2][0];
                    stream >> corners[2][1];

                    corners[1][0] = corners[2][0];corners[1][1] = corners[0][1];
                    corners[3][0] = corners[0][0];corners[3][1] = corners[2][1];

                    for(int k = 0 ; k < 4 ; k ++ ){
                        corners[k+4][0] = corners[k][0];
                        corners[k+4][1] = corners[k][1];
                    }
                    for(int  k = 0 ; k < 8 ; k++){
                        lead.shape.data[k] = QVector3D(corners[k][0],corners[k][1],corners[k][2]);
                    }
                    break;
                case(SHAPE_CONVEX_QUAD_XY):
                    break;
                case(SHAPE_RANGE_3D):
                    double x,y,z;

                    stream >> x;
                    stream >> y;
                    stream >> z;

                    QVector3D base = QVector3D(x,y,z);
                    stream >> x;
                    stream >> y;
                    stream >> z;
                    QVector3D dir  = QVector3D(x,y,z);
//                    qDebug() << base << dir;
                    // find max component of dir
                    int cpnt = 1 ;
                    if( dir.x() > dir.y() ){
                        if( dir.x() > dir.z() ) cpnt = 1;//x
                        else cpnt = 3; //z
                    }else if ( dir.y() > dir.z() ) cpnt = 2;
                    else cpnt = 3;
                    QVector3D someVec  = QVector3D( cpnt==3?1.0:0.0, cpnt==1?1.0:0.0 , cpnt==2?1.0:0.0 );
                    QVector3D tangent  = QVector3D::crossProduct(someVec,dir.normalized()).normalized();
                    QVector3D btangent = QVector3D::crossProduct(tangent,dir.normalized()).normalized();
//                    qDebug() << "someVec"  << someVec;
//                    qDebug() << "tangent"  << tangent;
//                    qDebug() << "btangent" << btangent;
                    double scale = (atoms_stats.max_corner -  atoms_stats.min_corner).length();
                    lead.shape.data[0] =  tangent*scale;
                    lead.shape.data[1] =  btangent*scale;
                    lead.shape.data[2] =  - tangent*scale;
                    lead.shape.data[3] =  - btangent*scale;

                    lead.shape.data[4] =   tangent*scale ;
                    lead.shape.data[5] =   btangent*scale;
                    lead.shape.data[6] =  - tangent*scale;
                    lead.shape.data[7] =  - btangent*scale;
//                    for(int  k = 0 ; k < 8 ; k++){
//                        if(lead.shape.data[k].x() < atoms_stats.min_corner.x()) lead.shape.data[k].setX(atoms_stats.min_corner.x());
//                        if(lead.shape.data[k].y() < atoms_stats.min_corner.y()) lead.shape.data[k].setY(atoms_stats.min_corner.y());
//                        if(lead.shape.data[k].z() < atoms_stats.min_corner.z()) lead.shape.data[k].setZ(atoms_stats.min_corner.z());

//                        if(lead.shape.data[k].x() > atoms_stats.max_corner.x()) lead.shape.data[k].setX(atoms_stats.max_corner.x());
//                        if(lead.shape.data[k].y() > atoms_stats.max_corner.y()) lead.shape.data[k].setY(atoms_stats.max_corner.y());
//                        if(lead.shape.data[k].z() > atoms_stats.max_corner.z()) lead.shape.data[k].setZ(atoms_stats.max_corner.z());
//                    }

                    lead.shape.data[0] += base;
                    lead.shape.data[1] += base;
                    lead.shape.data[2] += base;
                    lead.shape.data[3] += base;

                    lead.shape.data[4] += base + dir;
                    lead.shape.data[5] += base + dir;
                    lead.shape.data[6] += base + dir;
                    lead.shape.data[7] += base + dir;

                    break;
            }


        }
        // Next component

        // Check if the child tag name is COMPONENT
        if (Component.tagName()=="atoms")
        {

            // Get the first child of the component
            QDomElement Child=Component.firstChild().toElement();
            int id;
            // Read each child of the component node
            while (!Child.isNull())
            {
                // Read Name and value
                if (Child.tagName()=="d")   id=Child.firstChild().toText().data().toInt()-1;
                // Next child
                Child = Child.nextSibling().toElement();
                lead.atoms.push_back(id);
            }
        }
        // Check if the child tag name is COMPONENT
        if (Component.tagName()=="next_atoms")
        {

            // Get the first child of the component
            QDomElement Child=Component.firstChild().toElement();
            int id;
            // Read each child of the component node
            while (!Child.isNull())
            {
                // Read Name and value
                if (Child.tagName()=="d")   id=Child.firstChild().toText().data().toInt()-1;
                // Next child
                Child = Child.nextSibling().toElement();
                lead.nex_atoms.push_back(id);
            }
        }
        // Check if the child tag name is COMPONENT
        if (Component.tagName()=="lead_coupling")
        {

            // Get the first child of the component
            QDomElement Child=Component.firstChild().toElement();
            QString line;
            // Read each child of the component node
            while (!Child.isNull())
            {
                // Read Name and value
                if (Child.tagName()=="d")   line=Child.firstChild().toText().data();
                // Next child
                Child = Child.nextSibling().toElement();
                std::stringstream stream(line.toStdString());

                AtomConnection cnt;

                int a ;
                stream >> a;
                cnt.atomA = a-1;
                stream >> a;
                cnt.atomB = a-1;
                stream >> a;
                cnt.spinA = a;
                stream >> a;
                cnt.spinB = a;
                lead.cnts.push_back(cnt);
            }
        }
        // Check if the child tag name is COMPONENT
        if (Component.tagName()=="inner_coupling")
        {

            // Get the first child of the component
            QDomElement Child=Component.firstChild().toElement();
            QString line;
            // Read each child of the component node
            while (!Child.isNull())
            {
                // Read Name and value
                if (Child.tagName()=="d")   line=Child.firstChild().toText().data();
                // Next child
                Child = Child.nextSibling().toElement();
                std::stringstream stream(line.toStdString());

                AtomConnection cnt;

                int a ;
                stream >> a;
                cnt.atomA = a-1;
                stream >> a;
                cnt.atomB = a-1;
                stream >> a;
                cnt.spinA = a;
                stream >> a;
                cnt.spinB = a;
                lead.inner_cnts.push_back(cnt);
            }
        }
        // Next component
        Component = Component.nextSibling().toElement();
    }

    lead.visible = true;
    leads.push_back(lead);


}

void DataReader::precalculate_data(){
    QVector3D dims = atoms_stats.max_corner - atoms_stats.min_corner;

    QVector3D inv_dims = QVector3D(1/dims.x(),1/dims.y(),1/dims.z());

    if(dims.x() < 1.0e-10) inv_dims.setX(0.0);
    if(dims.y() < 1.0e-10) inv_dims.setY(0.0);
    if(dims.z() < 1.0e-10) inv_dims.setZ(0.0);

    double scale = 1.0;
    if(dims.x() > dims.y()){
        if(dims.x() > dims.z()) scale = inv_dims.x();
        else scale = inv_dims.z();
    }else if(dims.y() > dims.z()) scale = inv_dims.y();
    else scale = inv_dims.z();

    atoms_stats.scale = scale;
    atoms_stats.atom_radius = atoms_stats.ave_dist;


    p_atoms.clear();
    for(unsigned int i = 0 ; i < atoms.size() ; i++){
        Atom atom;
        atom     = atoms[i];
        atom.pos = (atoms[i].pos - atoms_stats.mass_center)*scale;
        p_atoms.push_back(atom);
    }
    p_connections.clear();
    for(unsigned int i = 0 ; i < connections.size() ; i++){
        AtomConnection& cnt = connections[i];
        if((cnt.spinA == atoms_stats.filter_spin_A && cnt.spinB == atoms_stats.filter_spin_B) ){
            p_connections.push_back(cnt);
        }
    }
    p_leads.clear();
    for(unsigned int l = 0 ; l < leads.size() ; l++){
        Lead& lead = leads[l];
        Lead p_lead;

        p_lead = lead;
        p_lead.cnts.clear();
        p_lead.inner_cnts.clear();

        for(int i = 0;i<12;i++){
            p_lead.shape.data[i] = (lead.shape.data[i] - atoms_stats.mass_center)*scale;
        }


        for(unsigned int i = 0 ; i < lead.cnts.size() ; i++){
            AtomConnection& cnt = lead.cnts[i];
            if((cnt.spinA == atoms_stats.filter_spin_A && cnt.spinB == atoms_stats.filter_spin_B) ){
                p_lead.cnts.push_back(cnt);
            }
        }

        for(unsigned int i = 0 ; i < lead.inner_cnts.size() ; i++){
            AtomConnection& cnt = lead.inner_cnts[i];
            if((cnt.spinA == atoms_stats.filter_spin_A && cnt.spinB == atoms_stats.filter_spin_B) ){
                p_lead.inner_cnts.push_back(cnt);
            }
        }
        p_leads.push_back(p_lead);
    }

}
