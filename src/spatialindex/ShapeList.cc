//
// Created by Chuang on 2019/6/13.
//
#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

ShapeList* ShapeList::clone() {
    throw Tools::NotSupportedException("clone");
}

ShapeList::ShapeList(const SpatialIndex::ShapeList &in) {
    m_MBRList=in.m_MBRList;
    m_MBCList=in.m_MBCList;
}
uint32_t ShapeList::getByteArraySize() const {
    int size= 3*sizeof(uint32_t);
    for(auto s:m_MBRList) size+=s->getByteArraySize();
    for(auto s:m_MBCList) size+=s->getByteArraySize();
    return size;
}

void ShapeList::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;
    uint8_t* tmpb;
    uint32_t tmplen;
    memcpy(ptr, &m_dimension, sizeof(uint32_t));
    memcpy(ptr, &m_datatype, sizeof(uint32_t));
    id_type size=(m_datatype==SpatialIndex::LeafBoundByMBR)?m_MBRList.size():m_MBCList.size();
    memcpy(ptr, &size, sizeof(id_type));
    ptr += sizeof(id_type);
    if(m_datatype==LeafBoundByMBR) {
        for (int i = 0; i < size; i++) {
            m_MBRList[i]->storeToByteArray(&tmpb, tmplen);
            memcpy(ptr, tmpb, tmplen);
            if (i != size - 1)
                ptr += tmplen;
        }
    }else if(m_datatype==LeafBoundByMBC){
        for (int i = 0; i < size; i++) {
            m_MBCList[i]->storeToByteArray(&tmpb, tmplen);
            memcpy(ptr, tmpb, tmplen);
            if (i != size - 1)
                ptr += tmplen;
        }
    }
}
void ShapeList::loadFromByteArray(const uint8_t *ptr) {
    memcpy(&m_dimension, ptr, sizeof(uint32_t));
    ptr+= sizeof(uint32_t);
    memcpy(&m_datatype, ptr, sizeof(uint32_t));
    ptr+= sizeof(uint32_t);
    id_type size;
    memcpy(&size, ptr, sizeof(id_type));
    ptr += sizeof(id_type);
    if(m_datatype==LeafBoundByMBR){
        m_MBRList.resize(size);
        uint32_t tmplen;
        for(int i=0;i<size;i++){
            m_MBRList[i]->loadFromByteArray(ptr);
            if(i!=size-1)
                ptr+=m_MBRList[i]->getByteArraySize();
        }
    }else if(m_datatype==LeafBoundByMBC){
        m_MBCList.resize(size);
        uint32_t tmplen;
        for(int i=0;i<size;i++){
            m_MBCList[i]->loadFromByteArray(ptr);
            if(i!=size-1)
                ptr+=m_MBCList[i]->getByteArraySize();
        }
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

void ShapeList::insert(SpatialIndex::IShape *shape) {
    if(m_datatype==-1) {
        Region *pmbr = dynamic_cast< Region *>(shape);
        if (pmbr != nullptr) {
            m_datatype = SpatialIndex::LeafBoundByMBR;
        }
        MBC *pmbc = dynamic_cast<MBC *>(shape);
        if (pmbc != nullptr) {
            m_datatype = SpatialIndex::LeafBoundByMBC;
        }
    }
    else if(m_datatype==SpatialIndex::LeafBoundByMBR){
        Region *pmbr = dynamic_cast< Region *>(shape);
        m_MBRList.push_back(pmbr);
    }
    else if(m_datatype==SpatialIndex::LeafBoundByMBC){
        MBC *pmbc = dynamic_cast<MBC *>(shape);
        m_MBCList.push_back(pmbc);
    }
}

void ShapeList::insert(SpatialIndex::Region *shape) {
    m_MBRList.push_back(shape);
}

void ShapeList::insert(SpatialIndex::MBC *shape) {
    m_MBCList.push_back(shape);
}