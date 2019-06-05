//
// Created by Chuang on 2019/5/16.
//

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include "spatialindex/tools/TimeDivision.h"

using std::string;
using std::vector;
//time division related


double naivetime(string l){
    int h;
    int m;
    int s;
    if(l.size()==8){
        h = stringToNum<int>(l.substr(0,2));
        m = stringToNum<int>(l.substr(3,5));
        s = stringToNum<int>(l.substr(6,8));
    } else{
        h = stringToNum<int>(l.substr(0,1));
        m = stringToNum<int>(l.substr(2,4));
        s = stringToNum<int>(l.substr(5,7));
    }
    return 3600*h+60*m+s;
}
int getPeriod(double time){
    if(time>getMaxPeriod()*PeriodLen)
        std::cerr<<"WARNING: calling get Period at "<<time<<" which is not supported. Maybe you should change the PeriodLen"<<std::endl;
    int pd= int(std::floor(time/PeriodLen));
    if(pd>=getMaxPeriod()) {
        pd=getMaxPeriod()-1;
    }
    return pd;
}
int getMaxPeriod(){//number of period, not id
    return 1;
}
double getPeriodStart(double time){
    int pd=getPeriod(time);
    return pd*PeriodLen;
}
double getPeriodEnd(double time){
    int pd=getPeriod(time);
    return (pd+1)*PeriodLen-0.00001;
}
