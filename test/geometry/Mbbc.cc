//
// Created by chuang on 4/2/19.
//

#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;
using namespace std;

int main(){
    double sLow[2]={0,0};
    double sHigh[2]={1,1};
    double eLow[2]={0,0};
    double eHigh[2]={1,1};
    double vLow[2]={1,1};
    double vHigh[2]={1,1};
    double pLow[2]={-1,2};
    double pHigh[2]={-1,2};
    double oLow[2]={1.3,1.3};
    double oHigh[2]={4,4};
    Mbbc m=*new Mbbc(*new Region(sLow,sHigh,2),*new Region(eLow,eHigh,2),
            *new Region(vLow,vHigh,2),*new Region(pLow,pHigh,2),0.0,1.0);
    TimeRegion t=*new TimeRegion(oLow,oHigh,0.5,0.5,2);
    cout<<m.intersectsTimeRegion(t);
    return 0;
}