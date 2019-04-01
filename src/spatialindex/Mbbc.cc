//
// Created by chuang on 4/1/19.
//

#include <cstring>
#include <cmath>
#include <limits>

#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

Mbbc::Mbbc(const SpatialIndex::Region &smbr, const SpatialIndex::Region &embr, const SpatialIndex::Region &vbr,
           const SpatialIndex::Region &pmbr, double tStart, double tEnd) {
    m_smbr=smbr;
    m_embr=embr;
    m_vmbr=vbr;
    m_pmbr=pmbr;
    m_startTime=tStart;
    m_endTime=tEnd;
}
Mbbc::Mbbc(const SpatialIndex::Mbbc &in) {
    m_smbr=in.m_smbr;
    m_embr=in.m_embr;
    m_vmbr=in.m_vmbr;
    m_pmbr=in.m_pmbr;
    m_startTime=in.m_startTime;
    m_endTime=in.m_endTime;
}

Mbbc& Mbbc::operator=(const Mbbc& r)
{
    if(this != &r)
    {
        m_smbr=r.m_smbr;
        m_embr=r.m_embr;
        m_vmbr=r.m_vmbr;
        m_pmbr=r.m_pmbr;
        m_startTime=r.m_startTime;
        m_endTime=r.m_endTime;
    }

    return *this;
}

bool Mbbc::operator==(const SpatialIndex::Mbbc &r) const {
    if (m_startTime < r.m_startTime - std::numeric_limits<double>::epsilon() ||
        m_startTime > r.m_startTime + std::numeric_limits<double>::epsilon() ||
        m_endTime < r.m_endTime - std::numeric_limits<double>::epsilon() ||
        m_endTime > r.m_endTime + std::numeric_limits<double>::epsilon())
        return false;
    if (!(m_smbr==r.m_smbr)||!(m_embr==r.m_embr)||!(m_vmbr==r.m_vmbr)||!(m_pmbr==r.m_pmbr))
        return false;
    return true;
}

Mbbc* Mbbc::clone() {
    return new Mbbc(*this);
}

