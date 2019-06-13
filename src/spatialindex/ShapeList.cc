//
// Created by Chuang on 2019/6/13.
//
#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

ShapeList* ShapeList::clone() {
    throw Tools::NotSupportedException("clone");
}

ShapeList::ShapeList(const SpatialIndex::ShapeList &in) {
    m_ShapeList=in.m_ShapeList;
}
uint32_t ShapeList::getByteArraySize() const {
    int size= 2*sizeof(uint32_t);
    for(auto s:m_ShapeList) size+=s->getByteArraySize();
    return size;
}

void ShapeList::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;
    uint8_t* tmpb;
    uint32_t tmplen;
    memcpy(ptr, &m_dimension, sizeof(uint32_t));
    id_type size=m_ShapeList.size();
    memcpy(ptr, &size, sizeof(id_type));
    ptr += sizeof(id_type);
    for(int i=0;i<size;i++){
        m_ShapeList[i]->storeToByteArray(&tmpb,tmplen);
        memcpy(ptr, tmpb, tmplen);
        if(i!=size-1)
            ptr += tmplen;
    }
}
void ShapeList::loadFromByteArray(const uint8_t *ptr) {
    memcpy(&m_dimension, ptr, sizeof(uint32_t));
    id_type size;
    memcpy(&size, ptr, sizeof(id_type));
    ptr += sizeof(id_type);
    m_ShapeList.resize(size);
    uint32_t tmplen;
    for(int i=0;i<size;i++){
        m_ShapeList[i]->loadFromByteArray(ptr);
        if(i!=size-1)
            ptr+=m_ShapeList[i]->getByteArraySize();
    }
}

//
// IEvolvingShape interface
//
void ShapeList::getVMBR(Region& out) const{
    throw Tools::NotSupportedException("clone");
}
void ShapeList::getMBRAtTime(double t, Region& out) const{
    throw Tools::NotSupportedException("clone");
}


//
// IShape interface
//
bool ShapeList::intersectsShape(const IShape& in) const{
    throw Tools::NotSupportedException("clone");
}
bool ShapeList::containsShape(const IShape& in) const{
    throw Tools::NotSupportedException("clone");
}
bool ShapeList::touchesShape(const IShape& in) const{
    throw Tools::NotSupportedException("clone");
}
void ShapeList::getCenter(Point& out) const{
    throw Tools::NotSupportedException("clone");
}
uint32_t ShapeList::getDimension() const{
    throw Tools::NotSupportedException("clone");
}
void ShapeList::getMBR(Region& out) const{
    throw Tools::NotSupportedException("clone");
}
void ShapeList::getTimeMBR(SpatialIndex::TimeRegion &out) const {
    throw Tools::NotSupportedException("clone");
}
double ShapeList::getArea() const{
    throw Tools::NotSupportedException("clone");
}
double ShapeList::getMinimumDistance(const IShape& in) const{
    throw Tools::NotSupportedException("clone");
}