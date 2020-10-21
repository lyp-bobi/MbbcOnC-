//
// Created by Chuang on 2020/10/18.
//

#ifndef SPATIALINDEX_GEOM_H
#define SPATIALINDEX_GEOM_H

#define precise double
#include <vector>
using namespace std;

class xPoint{
    precise x,y,t;
};
class xMBR{
    precise xmin,xmax,ymin,ymax,zmin,zmax;
};
class xMBC{
    xPoint ps,pe;
    precise rd,rv;
};
class xCylinder{
    xPoint c;
    precise r;
};
class xTraj{
    vector<xPoint> pts;
    bool fakehead, faketail;
};



#endif //SPATIALINDEX_GEOM_H
