////
//// Created by Chuang on 2019/6/12.
////
//
#include <spatialindex/SpatialIndex.h>
#include "xRTree.h"
#include "BulkLoader.h"
#include "storagemanager/xStore.h"
using namespace SpatialIndex::xRTreeNsp;


SpatialIndex::xRTreeNsp::xRTree * SpatialIndex::xRTreeNsp::createNewxRTree(IStorageManager *store, long indexfan, long leaffan, double bt) {
    Tools::Variant var;
    Tools::PropertySet ps;

    var.m_varType = Tools::VT_DOUBLE;
    var.m_val.dblVal = 0.99;
    ps.setProperty("FillFactor", var);

    var.m_varType = Tools::VT_DOUBLE;
    var.m_val.dblVal = bt;
    ps.setProperty("Bt", var);

    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = indexfan;
    ps.setProperty("IndexCapacity", var);

    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = leaffan;
    ps.setProperty("LeafCapacity", var);

    var.m_varType = Tools::VT_ULONG;
    var.m_val.ulVal = 3;
    ps.setProperty("Dimension", var);

    var.m_varType = Tools::VT_LONG;
    var.m_val.lVal = SpatialIndex::xRTreeNsp::RV_RSTAR;
    ps.setProperty("TreeVariant", var);

    SpatialIndex::xRTreeNsp::xRTree *ret = new SpatialIndex::xRTreeNsp::xRTree(*store, ps);

    ret->m_ts= static_cast<xStore *>(store);
//    var.m_varType = Tools::VT_LONGLONG;
//    var = ps.getProperty("IndexIdentifier");
//    indexIdentifier = var.m_val.llVal;

    return ret;
}

inline int idsize(){return 8;}
inline int nodeheadersize(){return 56;} //node id+mbr
inline int mbrsize(){return 48;}
inline int mbcsize(){
    if(bCompactMBC) return 56;
    else return 64;
}
inline int linesize(){return 48;}
inline int pointersize(){return 16;}
inline int entrysize(){return 32;}

xRTree * xRTreeNsp::buildMBRRTreeWP(IStorageManager *st,
                                    const CUTFUNC_PARA &f, double len, string add) {
    auto store=static_cast<xStore*>(st);
    auto stream = new xSBBStream(store,[&f,&len](auto &x){return f(x,len);});
    string name ="MBRWP"+std::to_string(len)+add;
    xRTree * r;
    tjstat->bt = len;
    if(store->m_property.contains(name)){
        Tools::Variant var;
        Tools::PropertySet ps;
        id_type id = store->m_property[name];
        var.m_varType = Tools::VT_LONGLONG;
        var.m_val.llVal = id;
        ps.setProperty("IndexIdentifier", var);
        r = new xRTree(*store, ps);
        r->m_bUsingMBR = true;
        r->m_ts=store;
        //std::cerr<<"load existing "<<name<<"\n";
    }
    else {
        std::cerr<<"start building "<<name<<"\n";
        int bindex = (PageSizeDefault - nodeheadersize()) / (idsize() + mbrsize()),
                bleaf = (PageSizeDefault - nodeheadersize()) / (idsize() + mbrsize() + pointersize() + entrysize());
        r = createNewxRTree(store, bindex, bleaf, len);
        r->m_bUsingMBR = true;
        store->m_property[name] = r->m_headerID;
        BulkLoader bl;
        bl.bulkLoadUsingSTR(r, *stream, bindex, bleaf, PageSizeDefault , 2000);
        std::cerr<<"built new "<<name<<"\n";
        st->flush();
    }
//    std::cerr<<r->m_headerID<<" "<<r->m_rootID<<endl;
    delete stream;
    return r;
}


xRTree * xRTreeNsp::buildMBCRTreeWP(IStorageManager *st,
                                    const CUTFUNC_PARA &f, double len, string add) {
    auto store=static_cast<xStore*>(st);

    auto stream = new xSBBStream(store,[&f,&len](auto &x){return f(x,len);});
    string name ="MBCWP"+std::to_string(len)+add;
    xRTree * r;
    tjstat->bt = len;
    if(store->m_property.contains(name)){
        Tools::Variant var;
        Tools::PropertySet ps;
        id_type id = store->m_property[name];
        var.m_varType = Tools::VT_LONGLONG;
        var.m_val.llVal = id;
        ps.setProperty("IndexIdentifier", var);
        r = new xRTree(*store, ps);
        r->m_bUsingMBC = true;
        r->m_ts=store;
//        std::cerr<<"load existing "<<name<<"\n";
    }
    else {
        std::cerr<<"start building "<<name<<"\n";
        int bindex = (PageSizeDefault - nodeheadersize()) / (idsize() + mbcsize()),
                bleaf = (PageSizeDefault - nodeheadersize()) / (idsize() + mbcsize() + pointersize() + entrysize());
        r = createNewxRTree(store, bindex, bleaf, len);
        r->m_bUsingMBC = true;
        store->m_property[name] = r->m_headerID;
        BulkLoader bl;
        bl.bulkLoadUsingSTR(r, *stream, bindex, bleaf, PageSizeDefault , 2000);
        std::cerr<<"built new "<<name<<"\n";
        st->flush();
    }
//    std::cerr<<r->m_headerID<<" "<<r->m_rootID<<endl;
    delete stream;
    return r;
}

xRTree * xRTreeNsp::buildTBTreeWP(IStorageManager *mng, string add) {
    auto store=static_cast<xStore*>(mng);
    auto stream = new xSBBStream(store, [](auto x){return xTrajectory::FP(x, 169);});
    string name ="TBWP"+add;
    xRTree * r;
    tjstat->bt = tjstat->tl*170;
    if(store->m_property.contains(name)){
        Tools::Variant var;
        Tools::PropertySet ps;
        id_type id = store->m_property[name];
        var.m_varType = Tools::VT_LONGLONG;
        var.m_val.llVal = id;
        ps.setProperty("IndexIdentifier", var);
        r = new xRTree(*mng,ps);
        r->m_bUsingMBR = true;
        r->m_ts=store;
//        std::cerr<<"load existing "<<name<<"\n";
    }
    else {
        std::cerr<<"start building "<<name<<"\n";
        int bindex = (PageSizeDefault - nodeheadersize()) / (idsize() + mbrsize()),
                bleaf = (PageSizeDefault - nodeheadersize()) / (idsize() + mbrsize() + pointersize() + entrysize());
        r = createNewxRTree(store, bindex, bleaf, tjstat->tl*170);
        r->m_bUsingMBR = true;
        store->m_property[name] = r->m_headerID;
        BulkLoader bl;
        bl.bulkLoadUsingSTR(r, *stream, bindex, bleaf, PageSizeDefault , 2000);
        std::cerr<<"built new "<<name<<"\n";
        mng->flush();
    }
    r->m_bUsingSBBD=r->m_bUsingSimp=false;
    delete stream;
    return r;
}

xRTree * xRTreeNsp::buildSTRTreeWP(IStorageManager *mng, string add) {
    auto store=static_cast<xStore*>(mng);
    auto stream = new xSBBStream(store, xTrajectory::EveryLine);
    string name ="STRWP"+add;
    xRTree * r;
    tjstat->bt = tjstat->tl;
    if(store->m_property.contains(name)){
        Tools::Variant var;
        Tools::PropertySet ps;
        id_type id = store->m_property[name];
        var.m_varType = Tools::VT_LONGLONG;
        var.m_val.llVal = id;
        ps.setProperty("IndexIdentifier", var);
        r = new xRTree(*mng,ps);
        r->m_bUsingMBL = true;
        r->m_ts=store;
//        std::cerr<<"load existing "<<name<<"\n";
    }
    else {
        std::cerr<<"start building "<<name<<"\n";
        int bindex = (PageSizeDefault - nodeheadersize()) / (idsize() + mbrsize()),
                bleaf = (PageSizeDefault - nodeheadersize()) / (idsize() + linesize() + pointersize() + entrysize());
        r = createNewxRTree(store, bindex, bleaf, tjstat->tl);
        r->m_bUsingMBL = true;
        store->m_property[name] = r->m_headerID;
        BulkLoader bl;
        bl.bulkLoadUsingSTR(r, *stream, bindex, bleaf, PageSizeDefault , 2000);
        std::cerr<<"built new "<<name<<"\n";
        mng->flush();
    }
    delete stream;
    r->m_bUsingSimp=false;
    return r;
}


xRTree * xRTreeNsp::buildMBRRTreeWoP(IStorageManager *st,
                                     const CUTFUNC_PARA &f, double len, string add) {
    auto store=static_cast<xStore*>(st);

    auto stream = new xSBBStream(store,[&f,&len](auto &x){return f(x,len);});
    string name ="MBRWoP"+std::to_string(len)+add;
    tjstat->bt = len;
    xRTree * r;
    if(store->m_property.contains(name)){
        Tools::Variant var;
        Tools::PropertySet ps;
        id_type id = store->m_property[name];
        var.m_varType = Tools::VT_LONGLONG;
        var.m_val.llVal = id;
        ps.setProperty("IndexIdentifier", var);
        r = new xRTree(*store, ps);
        r->m_bUsingMBR = true;
        r->m_bStoringLinks = false;
        r->m_ts=store;
        std::cerr<<"load existing "<<name<<"\n";
    }
    else {
        int bindex = (PageSizeDefault - nodeheadersize()) / (idsize() + mbrsize()),
                bleaf = (PageSizeDefault - nodeheadersize()) / (idsize() + mbrsize() + entrysize());
        r = createNewxRTree(store, bindex, bleaf, len);
        r->m_bUsingMBR = true;
        r->m_bStoringLinks = false;
        store->m_property[name] = r->m_headerID;
        BulkLoader bl;
        bl.bulkLoadUsingSTR(r, *stream, bindex, bleaf, PageSizeDefault , 2000);
        std::cerr<<"built new "<<name<<"\n";
        st->flush();
    }
//    std::cerr<<r->m_headerID<<" "<<r->m_rootID<<endl;
    delete stream;
    return r;
}



xRTree * xRTreeNsp::buildMBCRTreeWoP(IStorageManager *st,
                                     const CUTFUNC_PARA &f, double len, string add) {
    auto store=static_cast<xStore*>(st);

    auto stream = new xSBBStream(store,[&f,&len](auto &x){return f(x,len);});
    string name ="MBCWoP"+std::to_string(len) +add;
    xRTree * r;
    tjstat->bt = len;
    if(store->m_property.contains(name)){
        Tools::Variant var;
        Tools::PropertySet ps;
        id_type id = store->m_property[name];
        var.m_varType = Tools::VT_LONGLONG;
        var.m_val.llVal = id;
        ps.setProperty("IndexIdentifier", var);
        r = new xRTree(*store, ps);
        r->m_bUsingMBC = true;
        r->m_bStoringLinks = false;
        r->m_ts=store;
        std::cerr<<"load existing "<<name<<"\n";
    }
    else {
        int bindex = (PageSizeDefault - nodeheadersize()) / (idsize() + mbcsize()),
                bleaf = (PageSizeDefault - nodeheadersize()) / (idsize() + mbcsize() + entrysize());
        r = createNewxRTree(store, bindex, bleaf, len);
        r->m_bUsingMBC = true;
        r->m_bStoringLinks = false;
        store->m_property[name] = r->m_headerID;
        BulkLoader bl;
        bl.bulkLoadUsingSTR(r, *stream, bindex, bleaf, PageSizeDefault , 2000);
        std::cerr<<"built new "<<name<<"\n";
        st->flush();
    }
//    std::cerr<<r->m_headerID<<" "<<r->m_rootID<<endl;
    delete stream;
    return r;
}

xRTree * xRTreeNsp::buildTBTreeWoP(IStorageManager *mng) {
    auto store=static_cast<xStore*>(mng);
    auto stream = new xSBBStream(store, [](auto x){return xTrajectory::FP(x, 169);});
    string name ="TBWoP";
    xRTree * r;
    tjstat->bt = tjstat->tl*170;
    if(store->m_property.contains(name)){
        Tools::Variant var;
        Tools::PropertySet ps;
        id_type id = store->m_property[name];
        var.m_varType = Tools::VT_LONGLONG;
        var.m_val.llVal = id;
        ps.setProperty("IndexIdentifier", var);
        r = new xRTree(*mng,ps);
        r->m_bUsingMBR = true;
        r->m_bStoringLinks = false;
        r->m_ts=store;
        std::cerr<<"load existing "<<name<<"\n";
    }
    else {
        int bindex = (PageSizeDefault - nodeheadersize()) / (idsize() + mbrsize()),
                bleaf = (PageSizeDefault - nodeheadersize()) / (idsize() + mbrsize() + entrysize());
        r = createNewxRTree(store, bindex, bleaf, tjstat->tl * 170);
        r->m_bUsingMBR = true;
        r->m_bStoringLinks = false;
        store->m_property[name] = r->m_headerID;
        BulkLoader bl;
        bl.bulkLoadUsingSTR(r, *stream, bindex, bleaf, PageSizeDefault , 2000);
        std::cerr<<"built new "<<name<<"\n";
        mng->flush();
    }
    delete stream;
    return r;
}

xRTree * xRTreeNsp::buildSTRTreeWoP(IStorageManager *mng) {
    auto store=static_cast<xStore*>(mng);
    auto stream = new xSBBStream(store, xTrajectory::EveryLine);
    string name ="STRWoP";
    xRTree * r;
    tjstat->bt = tjstat->tl;
    if(store->m_property.contains(name)){
        Tools::Variant var;
        Tools::PropertySet ps;
        id_type id = store->m_property[name];
        var.m_varType = Tools::VT_LONGLONG;
        var.m_val.llVal = id;
        ps.setProperty("IndexIdentifier", var);
        r = new xRTree(*mng,ps);
        r->m_bUsingMBL = true;
        r->m_bStoringLinks = false;
        r->m_ts=store;
        std::cerr<<"load existing "<<name<<"\n";
    }
    else {
        int bindex = (PageSizeDefault - nodeheadersize()) / (idsize() + mbrsize()),
                bleaf = (PageSizeDefault - nodeheadersize()) / (idsize() + linesize() + entrysize());
        r = createNewxRTree(store, bindex, bleaf, tjstat->tl);
        r->m_bUsingMBL = true;
        r->m_bStoringLinks = false;
        store->m_property[name] = r->m_headerID;
        BulkLoader bl;
        bl.bulkLoadUsingSTR(r, *stream, bindex, bleaf, PageSizeDefault , 2000);
        std::cerr<<"built new "<<name<<"\n";
        mng->flush();
    }
    delete stream;
    return r;
}



xRTree * xRTreeNsp::loadTree(IStorageManager *mng, string name) {
    auto store=static_cast<xStore*>(mng);
    auto stream = new xSBBStream(store, xTrajectory::EveryLine);
    xRTree * r=NULL;
    if(store->m_property.contains(name)){
        Tools::Variant var;
        Tools::PropertySet ps;
        id_type id = store->m_property[name];
        var.m_varType = Tools::VT_LONGLONG;
        var.m_val.llVal = id;
        ps.setProperty("IndexIdentifier", var);
        r = new xRTree(*mng,ps);
        if(name.find("MBC")!=name.npos){
            r->m_bUsingMBC=true;
        }
        if(name.find("STR")!=name.npos){
            r->m_bUsingMBL=true;
        }
        if(name.find("WP")!=name.npos){
            r->m_bStoringLinks = true;
        }else{
            r->m_bStoringLinks = false;
        }
        r->m_ts=store;
        std::cerr<<"load existing "<<name<<"\n";
    }else{
        cerr<<"could not load the tree with given name"<<name<<endl;
    }
    delete stream;
    return r;
}

SBBForest * xRTreeNsp::buildSBBForest(IStorageManager *store, const function<queue<pair<pair<int, int>, xSBB> > (xTrajectory &, double)> &f, map<pair<double, double>, double> &lens,double slab) {
    SBBForest* res = new SBBForest();
    for(auto &i:lens){
        xRTree *r;
        if(i.second<=slab) {
            r = buildMBCRTreeWP(store, f, i.second);
        }else{
            r = buildMBRRTreeWP(store, f, i.second);
        }
        res->m_trees[i.first]=r;
    }
    return res;
}