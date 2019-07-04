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
    int size= 2*sizeof(uint32_t)+sizeof(id_type);
    for(const auto &s:m_MBRList) size+=s->getByteArraySize();
    for(const auto &s:m_MBCList) size+=s->getByteArraySize();
    return size;
}

void ShapeList::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;
    uint8_t* tmpb;
    uint32_t tmplen;
    memcpy(ptr, &m_dimension, sizeof(uint32_t));
    memcpy(ptr, &m_shapeType, sizeof(uint32_t));
    id_type size=(m_shapeType==SpatialIndex::LeafBoundByMBR)?m_MBRList.size():m_MBCList.size();
    memcpy(ptr, &size, sizeof(id_type));
    ptr += sizeof(id_type);
    if(m_shapeType==LeafBoundByMBR) {
        for (int i = 0; i < size; i++) {
            m_MBRList[i]->storeToByteArray(&tmpb, tmplen);
            memcpy(ptr, tmpb, tmplen);
            if (i != size - 1)
                ptr += tmplen;
        }
    }else if(m_shapeType==LeafBoundByMBC){
        for (int i = 0; i < size; i++) {
            m_MBCList[i]->storeToByteArray(&tmpb, tmplen);
            memcpy(ptr, tmpb, tmplen);
            if (i != size - 1)
                ptr += tmplen;
        }
    }
    else{
        std::cerr<<"bad shapeList\n";
    }
    std::cerr<<"len "<<len<<",expect"<<ptr-*data+tmplen<<"\n";
}
void ShapeList::loadFromByteArray(const uint8_t *ptr) {
    memcpy(&m_dimension, ptr, sizeof(uint32_t));
    ptr+= sizeof(uint32_t);
    memcpy(&m_shapeType, ptr, sizeof(uint32_t));
    ptr+= sizeof(uint32_t);
    id_type size;
    memcpy(&size, ptr, sizeof(id_type));
    ptr += sizeof(id_type);
    if(m_shapeType==LeafBoundByMBR){
        m_MBRList.resize(size);
        for(int i=0;i<size;i++){
            m_MBRList[i]->loadFromByteArray(ptr);
            if(i!=size-1)
                ptr+=m_MBRList[i]->getByteArraySize();
        }
    }else if(m_shapeType==LeafBoundByMBC){
        m_MBCList.resize(size);
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

//void ShapeList::insert(SpatialIndex::IShape *shape) {
//    if(m_shapeType==-1) {
//        Region *pmbr = dynamic_cast< Region *>(shape);
//        if (pmbr != nullptr) {
//            m_shapeType = SpatialIndex::LeafBoundByMBR;
//        }
//        MBC *pmbc = dynamic_cast<MBC *>(shape);
//        if (pmbc != nullptr) {
//            m_shapeType = SpatialIndex::LeafBoundByMBC;
//        }
//    }
//    else if(m_shapeType==SpatialIndex::LeafBoundByMBR){
//        Region *pmbr = dynamic_cast< Region *>(shape);
//        m_MBRList.emplace_back(pmbr);
//    }
//    else if(m_shapeType==SpatialIndex::LeafBoundByMBC){
//        MBC *pmbc = dynamic_cast<MBC *>(shape);
//        m_MBCList.emplace_back(pmbc);
//    }
//}

void ShapeList::insertbr(SpatialIndex::Region *shape) {
    m_MBRList.emplace_back(shape);
    m_shapeType=SpatialIndex::LeafBoundByMBR;
}

void ShapeList::insertbc(SpatialIndex::MBC *shape) {
    m_MBCList.emplace_back(shape);
    m_shapeType=SpatialIndex::LeafBoundByMBC;
}

std::ostream& SpatialIndex::operator<<(std::ostream& os, const ShapeList& r)
{
    if(r.m_shapeType==LeafBoundByMBR){
        os<<"Shape List have "<<r.m_MBRList.size()<<" MBRs.\n";
        for(auto br:r.m_MBRList){
            os<<*br;
        }
    }
    else if(r.m_shapeType==LeafBoundByMBC){
        os<<"Shape List have "<<r.m_MBCList.size()<<" MBCs.\n";
        for(auto bc:r.m_MBCList){
            os<<*bc;
        }
    }
    os<<"\n";
	return os;
}
