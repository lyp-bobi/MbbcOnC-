//
// Created by Chuang on 2021/1/27.
//

#include "testFuncs.h"

int main(){
    cout<<"specify a file"<<endl;
    string target;
    cin>>target;
    xStore x(target, testFileName(target), true);
    xRTree * r=NULL;
    cout<<"loaded store\n";
    cout<<x.m_property<<endl;
    char command;


    while(cin>>command){
        if(command=='r'){
            cout<<"give tree name\n";
            string para;
            cin>>para;
            if(r!=NULL) delete r;
            r=loadTree(&x, para);
            if(r!=NULL){
                cout<<"tree root"<<r->m_rootID<<endl;
            }
        }else if(command=='n'){
            cout<<"give node number:\n";
            id_type idp;
            cin>>idp;
            auto n=r->readNode(idp);
            cout<<n->toString()<<endl;
        }else if(command=='t'){
            cout<<"give traj id:\n";
            xTrajectory s;
            id_type idp;
            int ps,pe;
            cin>>idp>>ps>>pe;
            x.loadTraj(s,xStoreEntry(idp,ps,pe));
            cout<<s.toString()<<endl;
        }
    }

}