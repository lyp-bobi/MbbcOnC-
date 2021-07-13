//
// Created by Chuang on 2021/4/8.
//
#include <iostream>
using std::cerr;
using std::endl;
#include <sqlite3.h>
#include "../../include/storagemanager/tjsql.h"
#include <thread>


/* todo: free these stmts*/
thread_local sqlite3* db = NULL;
static string dbfile;
thread_local sqlite3_stmt *stmt_page_insert = NULL;
thread_local sqlite3_stmt *stmt_page_load = NULL;
thread_local sqlite3_stmt *stmt_traj_insert = NULL;
thread_local sqlite3_stmt *stmt_traj_load = NULL;

thread_local db_thread dbt;
static int64_t next_page = -1;
static int64_t next_traj = -1;

db_thread::db_thread() {
}

void db_thread::release() {
    if(db != NULL){
        sqlite3_finalize(stmt_page_insert);
        sqlite3_finalize(stmt_page_load);
        sqlite3_finalize(stmt_traj_insert);
        sqlite3_finalize(stmt_traj_load);
        sqlite3_close(db);
        db = NULL;
    }
}

db_thread::~db_thread() noexcept {
    if(db != NULL){
        sqlite3_finalize(stmt_page_insert);
        sqlite3_finalize(stmt_page_load);
        sqlite3_finalize(stmt_traj_insert);
        sqlite3_finalize(stmt_traj_load);
        sqlite3_close(db);
        db = NULL;
    }
}

bool conn_init(string path){
    int rc;
    dbfile = path;
    if(db==NULL){
        rc = sqlite3_open(path.c_str(), &db);
        if (rc != SQLITE_OK) {
            cerr << "db open failed: " << sqlite3_errmsg(db) << endl;
        }
        sqlite3_busy_timeout(db, 30000);
        db_create_table();
        string sql;
        sql ="INSERT OR REPLACE INTO page(PAGEID, DATA)"
                                       " VALUES(?1,?2)";
        rc = sqlite3_prepare(db,sql.c_str(),-1, &stmt_page_insert, NULL);
        if (rc != SQLITE_OK) {
            cerr << sqlite3_errmsg(db) << endl;
        }
        sql = "SELECT DATA FROM page WHERE PAGEID = ?1";
        rc = sqlite3_prepare(db,sql.c_str(),-1, &stmt_page_load, NULL);
        if (rc != SQLITE_OK) {
            cerr << sqlite3_errmsg(db) << endl;
        }
        sql ="INSERT INTO traj(TRAJID, PAGEID, LENGTH)"
                                " VALUES(?1,?2,?3)";
        rc = sqlite3_prepare(db,sql.c_str(),-1, &stmt_traj_insert, NULL);
        if (rc != SQLITE_OK) {
            cerr << sqlite3_errmsg(db) << endl;
        }
        sql = "SELECT PAGEID, LENGTH FROM traj WHERE TRAJID = ?1";
        rc = sqlite3_prepare(db,sql.c_str(),-1, &stmt_traj_load, NULL);
        if (rc != SQLITE_OK) {
            cerr << sqlite3_errmsg(db) << endl;
        }
    }
    return true;
}

bool db_create_table(){
    int rc;
    string sql =
            "create table if not exists page(PAGEID INTEGER PRIMARY KEY ASC, DATA BLOB)";
    rc = sqlite3_exec(db,sql.c_str(),NULL,NULL,NULL);
    if (rc != SQLITE_OK) {
        cerr << "create table " << sqlite3_errmsg(db) << endl;
    }
    sql =
            "create table if not exists traj(TRAJID INTEGER PRIMARY KEY ASC, PAGEID INT8, LENGTH INT4)";
    rc = sqlite3_exec(db,sql.c_str(),NULL,NULL,NULL);
    if (rc != SQLITE_OK) {
        cerr<<"create table " << sqlite3_errmsg(db) << endl;
    }
    return true;
}

bool db_insert_page(int64_t &pageid, uint32_t len, const uint8_t* const data){
    if(pageid == -1) {
        pageid = db_last_pageid() + 1;
    }
    int rc = SQLITE_ABORT;
    while(rc !=SQLITE_DONE) {
        rc = sqlite3_reset(stmt_page_insert);
        rc = sqlite3_bind_int64(stmt_page_insert, 1, pageid);
        rc = sqlite3_bind_blob(stmt_page_insert, 2, data, len, SQLITE_STATIC);
        rc = sqlite3_step(stmt_page_insert);
        if (rc != SQLITE_DONE) {
            cerr << "insert page: " << sqlite3_errmsg(db) << endl;
        }
    }
    next_page += 1;
    return true;
}

bool db_load_page(int64_t pageid, uint32_t& len, uint8_t** data){
    int rc;
    if(db == NULL) conn_init(dbfile);
    rc = sqlite3_reset(stmt_page_load);
    rc = sqlite3_bind_int64(stmt_page_load,1,pageid);
    rc = sqlite3_step(stmt_page_load);
    if(rc!= SQLITE_ROW) {
        cerr << "load page " << pageid << "failed";
        cerr << sqlite3_errmsg(db) << endl;
    }
    *data = (unsigned char*) sqlite3_column_blob(stmt_page_load,0);
    len = sqlite3_column_bytes(stmt_page_load, 0);
    return true;
}

int64_t db_last_pageid(void){
    if(db == NULL) conn_init(dbfile);
    if(next_page>=0)
        return next_page;
    int rc;
    sqlite3_stmt *stmt;
    string sql;
    sql ="SELECT max(PAGEID) FROM page";
    rc = sqlite3_prepare(db,sql.c_str(),-1, &stmt, NULL);
    rc = sqlite3_step(stmt);
    if(rc != SQLITE_ROW)
    {
        cerr<< "last page id: " << sqlite3_errmsg(db) << endl;
        return 0;
    }
    int64_t res = sqlite3_column_int64(stmt, 0);
    sqlite3_finalize(stmt);
    next_page = res;
    return res;
}

int64_t db_last_trajid(void){
    if(db == NULL) conn_init(dbfile);
    if(next_traj>=0)
        return next_traj;
    int rc;
    sqlite3_stmt *stmt;
    string sql;
    sql ="SELECT max(TRAJID) FROM traj";
    rc = sqlite3_prepare(db,sql.c_str(),-1, &stmt, NULL);
    rc = sqlite3_step(stmt);
    if(rc != SQLITE_ROW)
    {
        cerr<< "last traj id: " << sqlite3_errmsg(db) << endl;
        return 0;
    }
    int64_t res = sqlite3_column_int64(stmt, 0);
    sqlite3_finalize(stmt);
    next_traj = res;
    return res;
}

bool db_insert_traj(int64_t id, int64_t pageid, uint32_t npoint){
    if(db == NULL) conn_init(dbfile);
    int rc;
    rc = sqlite3_reset(stmt_traj_insert);
    rc = sqlite3_bind_int64(stmt_traj_insert,1,(int64_t)id);
    rc = sqlite3_bind_int64(stmt_traj_insert,2,(int64_t)pageid);
    rc = sqlite3_bind_int(stmt_traj_insert,3,(int32_t)npoint);
    rc = sqlite3_step(stmt_traj_insert);
    if (rc != SQLITE_DONE) {
        cerr <<"insert traj:" << sqlite3_errmsg(db) << endl;
    }
    return true;
}

xTrajEntry db_load_traj_entry(int64_t trajid){
    if(db == NULL) conn_init(dbfile);
    int rc;
    xTrajEntry res;
    rc = sqlite3_reset(stmt_traj_load);
    rc = sqlite3_bind_int64(stmt_traj_load,1,trajid);
    rc = sqlite3_step(stmt_traj_load);
    if (rc != SQLITE_ROW) {
        cerr <<"load traj entry:" << sqlite3_errmsg(db) << endl;
    }
    res.m_page = sqlite3_column_int64(stmt_traj_load,0);
    res.m_npoint = sqlite3_column_int(stmt_traj_load,1);
    return res;
}