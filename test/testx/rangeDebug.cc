//
// Created by Chuang on 2021/2/4.
//

#include "testFuncs.h"
#include "random"

int main(){
    try {
        string target = "tdfilter.txt";
        double len = 600;
        cerr<<"seglen: ";
        cerr<<endl;
        xStore x(target, testFileName(target), true);
        vector<xCylinder> queries;
        fillQuerySet(queries,x,0.1,0);
        xCylinder query(xPoint(116.327,40,6516),4,226516,226516,2);
        xTrajectory tj2;
        int sum = 0;
        for(auto &s:*(x.m_trajIdx)){
            x.loadTraj(tj2,xStoreEntry(s.first,0,100000));
            if(tj2.intersectsxCylinder(query)){
                sum +=1;
            }else{
                cerr<<s.first<<"\t"<<tj2<<endl;
                cerr<<tj2.intersectsxCylinder(query)<<endl;
            }
        }
        cerr<<sum<<endl;
        auto r = buildMBCRTreeWP(&x, xTrajectory::ISS, len);
        MyVisitor vis;
        r->intersectsWithQuery(query,vis);
        cerr<<vis.m_resultGet<<endl;
        cerr<<"mission complete.\n";
    }catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    return 0;
}

