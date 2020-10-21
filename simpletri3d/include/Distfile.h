//
// Created by Chuang on 2020/10/18.
//

#ifndef SPATIALINDEX_DISKFILE_H
#define SPATIALINDEX_DISKFILE_H
#include <stdint.h>
#include <vector>
#include <fstream>
#include <set>
#include <map>
#define id_type long long

class DiskStorageManager
{
public:
    DiskStorageManager();
    virtual ~DiskStorageManager();

    virtual void flush();
    virtual void loadByteArray(const id_type page, uint32_t& len, uint8_t** data);
    virtual void storeByteArray(id_type& page, const uint32_t len, const uint8_t* const data);
    virtual void deleteByteArray(const id_type page);

private:
    class Entry
    {
    public:
        uint32_t m_length;
        std::vector<id_type> m_pages;
    };

protected:
    std::fstream m_dataFile;
    std::fstream m_indexFile;
    uint32_t m_pageSize;
    id_type m_nextPage;
    std::set<id_type> m_emptyPages;
    std::map<id_type, Entry*> m_pageIndex;

    uint8_t* m_buffer;
public:
    double iotime=0;
    id_type nextPage() {return m_nextPage;}
}; // DiskStorageManager

#endif //SPATIALINDEX_DISKFILE_H
