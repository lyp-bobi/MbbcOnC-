//
// Created by Chuang on 2019/5/29.
//
#define _USE_MATH_DEFINES

#include <cstring>
#include <math.h>
#include <limits>
#include <algorithm>
#include <cfloat>

#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;
xSBB::xSBB(){}
xSBB::xSBB(const SpatialIndex::xSBB &in)
:hasbr(in.hasbr),hasbc(in.hasbc),hasbl(in.hasbl){
    if(hasbr){
        br=new xMBR(*in.br);
    }
    if(hasbc){
        bc=new xMBC(*in.bc);
    }
    if(hasbl){
        bl=new xLine(*in.bl);
    }
}
xSBB::xSBB(xMBR &r){
    hasbr=true;
    br=new xMBR(r);
}
xSBB::xSBB(xMBC &r){
    hasbc=true;
    bc=new xMBC(r);
}
xSBB::xSBB(xLine &r){
    hasbl=true;
    bl=new xLine(r);
}
xSBB::~xSBB(){
    if(hasbr){
        delete br;
    }
    if(hasbc){
        delete bc;
    }
    if(hasbl){
        delete bl;
    }
}


xSBB& xSBB::operator=(const xSBB& in)
{
    hasbr=in.hasbr;hasbc=in.hasbc;hasbl=in.hasbl;
    if(hasbr){
        br=new xMBR(*in.br);
    }
    if(hasbc){
        bc=new xMBC(*in.bc);
    }
    if(hasbl){
        bl=new xLine(*in.bl);
    }
    return *this;
}

bool xSBB::operator==(const SpatialIndex::xSBB &r) const {
    throw Tools::IllegalArgumentException(
            "SBB:== not supported"
    );
    return true;
}
//
// IObject interface
//
xSBB* xSBB::clone() {
    return new xSBB(*this);
}

void xSBB::loadbr(xMBR &r) {
    if(hasbr) delete br;
    hasbr=true;
    br = new xMBR(r);
}

void xSBB::loadbc(xMBC &r) {
    if(hasbc) delete bc;
    hasbc=true;
    bc = new xMBC(r);
}

void xSBB::loadbl(xLine &r) {
    if(hasbl) delete bl;
    hasbl=true;
    bl = new xLine(r);
}

double xSBB::tdist(const xPoint &p) const {
    if(hasbr){
        return br->getMinimumDistance(p);
    }
    if(hasbc){
        return bc->getMinimumDistance(p);
    }
    if(hasbl){
        return bl->getMinimumDistance(p);
    }
    return 0;
}

double xSBB::startTime() const {
    if(hasbr){
        return br->m_tmin;
    }
    if(hasbc){
        return bc->m_ps.m_t;
    }
    if(hasbl){
        return bl->m_ps.m_t;
    }
    return 0;
}

double xSBB::endTime() const {
    if(hasbr){
        return br->m_tmax;
    }
    if(hasbc){
        return bc->m_pe.m_t;
    }
    if(hasbl){
        return bl->m_pe.m_t;
    }
    return 0;
}


std::string xSBB::toString() const {
    std::string s ="";
    if(hasbr) s += "1 ";
    if(hasbc) s += "1 ";
    if(hasbl) s += "1 ";
    if(hasbr){
        s += std::to_string(br->m_xmin) + " " + std::to_string(br->m_ymin) + " "
                +std::to_string(br->m_tmin);
        s += std::to_string(br->m_xmax) + " " + std::to_string(br->m_ymax) + " "
             +std::to_string(br->m_tmax);
    }
    if(hasbc) {
        s += std::to_string(bc->m_ps.m_x) + " " + std::to_string(bc->m_ps.m_y) + " " + std::to_string(bc->m_ps.m_t) + " ";
        s += std::to_string(bc->m_pe.m_x) + " " + std::to_string(bc->m_pe.m_y) + " " + std::to_string(bc->m_pe.m_t) + " ";
        s += std::to_string(bc->m_rd) + " " + std::to_string(bc->m_rv);
    }
    if(hasbc) {
        s += std::to_string(bc->m_ps.m_x) + " " + std::to_string(bc->m_ps.m_y) + " " + std::to_string(bc->m_ps.m_t) + " ";
        s += std::to_string(bc->m_pe.m_x) + " " + std::to_string(bc->m_pe.m_y) + " " + std::to_string(bc->m_pe.m_t) + " ";
    }
    return s;
}
void xSBB::loadFromString(std::string s) {
    auto nums=split(s,' ');
    int cur=0;
    if(std::stod(nums[cur++])==1) hasbr =true;
    if(std::stod(nums[cur++])==1) hasbc =true;
    if(std::stod(nums[cur++])==1) hasbl =true;
    if(hasbr){
        if(br!= nullptr) delete br;
        br=new xMBR();
        br->m_xmin=std::stod(nums[cur++]);
        br->m_ymin=std::stod(nums[cur++]);
        br->m_tmin=std::stod(nums[cur++]);
        br->m_xmax=std::stod(nums[cur++]);
        br->m_ymax=std::stod(nums[cur++]);
        br->m_tmax=std::stod(nums[cur++]);
    }
    if(hasbc) {
        if(bc!= nullptr) delete bc;
        bc=new xMBC();
        bc->m_ps.m_x = std::stod(nums[cur++]);
        bc->m_ps.m_y = std::stod(nums[cur++]);
        bc->m_ps.m_t = std::stod(nums[cur++]);
        bc->m_pe.m_x = std::stod(nums[cur++]);
        bc->m_pe.m_y = std::stod(nums[cur++]);
        bc->m_pe.m_t = std::stod(nums[cur++]);
        bc->m_rd = std::stod(nums[cur++]);
        bc->m_rv = std::stod(nums[cur++]);
    }
    if(hasbl) {
        if(bl!= nullptr) delete bl;
        bl=new xLine();
        bl->m_ps.m_x = std::stod(nums[cur++]);
        bl->m_ps.m_y = std::stod(nums[cur++]);
        bl->m_ps.m_t = std::stod(nums[cur++]);
        bl->m_pe.m_x = std::stod(nums[cur++]);
        bl->m_pe.m_y = std::stod(nums[cur++]);
        bl->m_pe.m_t = std::stod(nums[cur++]);
    }
}