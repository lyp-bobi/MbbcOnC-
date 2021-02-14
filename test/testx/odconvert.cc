//
// Created by Chuang on 2021/2/10.
//
#include "testFuncs.h"

int main(){
    vector<pair<id_type, xTrajectory> > res;
    vector<string> files;
    struct dirent *ptr;
    DIR *dir;
    string PATH = "D://out2/";
    dir = opendir(PATH.c_str());
    while ((ptr = readdir(dir)) != NULL) {
        if (ptr->d_name[0] == '.')
            continue;
        //cout << ptr->d_name << endl;
        files.emplace_back(PATH + ptr->d_name);
    }
    id_type id = -1;
    for(int i=0;i<files.size();i+=5){
        auto j=loadGTToTrajs(files[i]);
        id = dumpToFile(j,to_string(i/5)+".out",-1,id);
        j= loadGTToTrajs(files[i+1]);
        id = dumpToFile_append(j,to_string(i/5)+".out",-1,id);
        j= loadGTToTrajs(files[i+2]);
        id = dumpToFile_append(j,to_string(i/5)+".out",-1,id);
        j= loadGTToTrajs(files[i+3]);
        id =dumpToFile_append(j,to_string(i/5)+".out",-1,id);
        j= loadGTToTrajs(files[i+4]);
        id =id=dumpToFile_append(j,to_string(i/5)+".out",-1,id);
    }
    return 0;
}