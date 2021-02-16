//
// Created by Chuang on 2021/2/17.
//

#include "testFuncs.h"
#include "random"

int main(int argc,char *argv[]){
    try {
        string target;
        if(argc==1) {
            target = "tdexpand.data";
        }else {
            target = "glexpand.data";
        }
        double seglens[] = {600,900,1200,1800,2700,3600};
        cerr<<"seglen: ";
        for(auto len:seglens){cerr<<len<<" ";}
        cerr<<endl;
        cerr<<"12,sbbdknn"<<target<<endl;
        xStore x(target, testFileName(target), true);
        for(double qt=300;qt<=5300;qt+=500) {
            cerr<<"qt is " << qt<<endl;
            vector<xTrajectory> queries;
            fillQuerySet(queries,x,qt);
            for (auto len:seglens) {
                MTQ q;
                q.prepareTrees(&x, [&len](auto x) {
                    xRTree* r= buildMBCRTreeWP(x, xTrajectory::OPTS, len);
                    r->m_bUsingSBBD=false;
                    return r;
                }
                    );
                q.appendQueries(queries);
                std::cerr << q.runQueries().toString();
            }
            for (auto len:seglens) {
                MTQ q;
                q.prepareTrees(&x, [&len](auto x) {
                    return buildMBCRTreeWP(x, xTrajectory::OPTS, len);
                });
                q.appendQueries(queries);
                std::cerr << q.runQueries().toString();
            }
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