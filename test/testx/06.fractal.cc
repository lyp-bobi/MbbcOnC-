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
        for(double rd = 4;rd>0.001;rd/=2) {
            cerr<<"rd is " << rd<<endl;
            vector<xCylinder> queries;
            fillQuerySet(queries,x,rd,0);
            MTQ q;
            q.prepareTrees(&x, [&len](auto x) { return buildMBCRTreeWP(x, xTrajectory::OPTS, len); });
            q.appendQueries(queries);
            std::cerr << q.runQueries().toString();
        }
        cerr<<"mission complete.\n";
    }catch (Tools::Exception &e) {
        cerr << "******ERROR******" << endl;
        std::string s = e.what();
        cerr << s << endl;
        return -1;
    }
    return 0;
}