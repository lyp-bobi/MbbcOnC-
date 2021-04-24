//
// Created by Chuang on 2021/4/10.
//

#include <iostream>
#include <storagemanager/tjsql.h>

int main(){
    conn_init("sqldb");
    db_insert_traj(1,1,10);
    uint64_t p=50;
    uint32_t n=20;
    std::cerr<<p<<"   "<<n;
    return 0;
}