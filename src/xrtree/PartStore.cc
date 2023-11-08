//
// Created by Chuang on 2021/2/3.
//

#include "xRTree.h"

using namespace xRTreeNsp;


void PartsStore::Parts::insert(xSBB &r, id_type prev, id_type next, xStoreEntry &entry) {
    is_modified = true;
    m_UCsbbs.emplace_back(r);
    m_loadedTime += min(r.m_endTime, m_ps->m_simpquery.m_endTime()) -
                    max(r.m_startTime, m_ps->m_simpquery.m_startTime());
    if(r.m_endTime > m_maxtime){
        m_maxtime=min(r.m_endTime, m_ps->m_simpquery.m_endTime());
        if(next==-1|| r.m_endTime > m_ps->m_simpquery.m_endTime()) {
            m_hasNext=false;
            m_loadedTime+= m_ps->m_simpquery.m_endTime() - m_maxtime;
        }
    }
    if(r.m_startTime < m_mintime){
        m_mintime=max(r.m_startTime, m_ps->m_simpquery.m_startTime());
        if(prev==-1|| r.m_startTime < m_ps->m_simpquery.m_startTime()) {
            m_hasPrev=false;
            m_loadedTime+=m_mintime-m_ps->m_simpquery.m_startTime();
        }
    }
    if(prev>=0&&m_ps->m_loadedLeaf.find(prev)==m_ps->m_loadedLeaf.end()) m_missingLeaf.insert(prev);
    if(next>=0&&m_ps->m_loadedLeaf.find(next)==m_ps->m_loadedLeaf.end()) m_missingLeaf.insert(next);
    m_ses[r.m_startTime]=entry;
}

bool PartsStore::insert(id_type id, xSBB &b, id_type leafid, id_type prev, id_type next, xStoreEntry &entry) {
    bool res=false;
    auto it = m_parts.find(id);
    if(it==m_parts.end()){
        m_parts.emplace(make_pair(id,Parts(this)));
        it = m_parts.find(id);
        res = true;
    }
    it->second.insert(b,prev,next,entry);
    it->second.m_missingLeaf.erase(leafid);
    return res;
}

DISTE PartsStore::updateValue(id_type id,bool Inqueue) {
    Parts *parts = &m_parts[id];
//    if(!parts->is_modified) return parts->m_calcMin;
    std::pair<double, double> timeInterval;
    double timesum=0;
//    for(auto &b:parts->m_UCsbbs){
//        timesum +=b.m_endTime-b.m_startTime;
//    }
//    if(!parts->m_missingLeaf.empty()
//        &&timesum*5<m_simpquery.m_endTime() - m_simpquery.m_startTime())
//    {
//        return parts->m_calcMin;
//    }

    DISTE res;
    int type = 2;
    if (parts->m_missingLeaf.empty()) type = 3;

    if(current_distance != RMDTW && current_distance != IED)
    {
        return DISTE(0);
    }

    if(type == 3) {
        if (current_distance == RMDTW)
        {
            if (parts->m_UCsbbs.size() > 0) {
                parts->m_sbbs.insert(parts->m_sbbs.end(),
                                     parts->m_UCsbbs.begin(),
                                     parts->m_UCsbbs.end());
                parts->m_UCsbbs.clear();
            }
            if (parts->m_completelb.infer) {
                std::sort(parts->m_sbbs.begin(), parts->m_sbbs.end(),[](const xSBB &a, const xSBB &b)
                {
                    return a.m_startTime < b.m_startTime;
                });
                vector<pair<xPoint, double>> cross;
                int cur1 = 0, cur2 = 0;
                for(cur1 = 0; cur1 < m_simpquery.m_points.size();cur1++)
                {
                    while(parts->m_sbbs[cur2].m_endTime < m_simpquery.m_points[cur1].m_t
                        && cur2 != parts->m_sbbs.size() - 1)
                        cur2++;
                    cross.emplace_back(parts->m_sbbs[cur2].crossSec(m_simpquery.m_points[cur1].m_t));
                }
                parts->m_completelb = m_simpquery.getRMDTW(cross);
            }
            res = parts->m_completelb;
        }
    } else {
        for(auto &b:parts->m_UCsbbs){
            parts->putSBB(b);
        }
        parts->m_UCsbbs.clear();
        double slab = (m_simpquery.m_endTime() - m_simpquery.m_startTime()) /
                      (m_simpquery.m_points.size() - 1);

        for (auto it = parts->m_line.begin(); it != parts->m_line.end(); it++) {
            if (!it->d.infer) {
                res += it->d;
            } else {
                if (!m_nodespq.empty()) {
                    if (current_distance == RMDTW) {
                        res += DISTE(
                                m_nodespq.top()->m_dist.opt /
                                m_simpquery.m_points.size()
                                * ((int) ((it->te - it->ts) / slab)),
                                1e300, true);
                    } else if (current_distance == IED) {
                        if (!m_nodespq.empty())
                            res += DISTE(
                                    min(it->d.opt, (it->te - it->ts) *
                                                   m_nodespq.top()->m_dist.opt /
                                                   (m_simpquery.m_endTime() -
                                                    m_simpquery.m_startTime())),
                                    1e300, true);
                    }
                }
            }
        }
    }
    parts->m_calcMin = res;

    res.opt -= m_error;
    res.pes += m_error;
    m_pes.insert(id, res.pes);
    if(res.opt>m_pes.threshold()){
        m_except.insert(id);
    }
    if(Inqueue)
        m_mpq.updateValue(m_handlers[id], id, res, type);
    parts->is_modified=false;
    return res;
}


void PartsStore::loadLeaf(const Node &n, double dist) {
//                    std::cerr<<"load leaf"<<n.m_nodeMBR<<"\n";
//                    std::cerr<<"leaf dist"<<m_query.getNodeMinimumDistance(n.m_nodeMBR,100)/(m_query.m_endTime()-m_query.m_startTime())<<"\n";
//                    std::cerr<<"load leaf"<<n.m_identifier<<"\n";
    m_loadedLeaf.insert(n.m_identifier);
    std::set<id_type> nid;
    for(int i=0;i<n.m_children;i++){
        id_type trajid=n.m_se[i].m_id;
        xStoreEntry entry= n.m_se[i];
        double bts= n.m_ptrxSBB[i]->m_startTime,bte= n.m_ptrxSBB[i]->m_endTime;
        if(bts >= m_simpquery.m_endTime() ||
           bte <= m_simpquery.m_startTime() || m_except.find(trajid) != m_except.end()){}
        else {
            if(insert(trajid, *n.m_ptrxSBB[i],n.m_identifier,
                      (m_simpquery.m_startTime() < bts) ? n.m_prevNode[i] : -1
                    , (m_simpquery.m_endTime() > bte) ? n.m_nextNode[i] : -1,
                   entry)){
                nid.insert(trajid);
            }
        }
    }
    for(const auto &rid:nid){
        DISTE pd=updateValue(rid,false);
        auto handle = m_mpq.push(new NNEntry(rid, pd, 2+(pd.infer?0:1)));
        m_handlers[rid] = handle;
    }
}

NNEntry* PartsStore::top() {
    id_type lastid = -1;
    if(m_mpq.empty()&&m_nodespq.empty()) return NULL;
    while(true){
        if(m_mpq.empty()||(!m_nodespq.empty()&&m_nodespq.top()->m_dist < m_mpq.top()->m_dist)){
            return m_nodespq.top();
        }
        if(m_mpq.top()->m_type>=3)
            return m_mpq.top();
        if(lastid!=m_mpq.top()->m_id) {
            if(m_except.find(m_mpq.top()->m_id)!=m_except.end()){
                NNEntry* p = m_mpq.top();
                m_mpq.pop();
                delete p;
                lastid = -1;
            }else {
                lastid = m_mpq.top()->m_id;
                updateValue(lastid);
                m_mpq.updateOrder(m_handlers[lastid]);
            }
        } else{
            return m_mpq.top();
        }
    }
}


NNEntry* PartsStore::pop(int type) {
    if(type <2)
        return m_nodespq.pop();
    else
        return m_mpq.pop();
}

void PartsStore::push(NNEntry *e) {
    if(e->m_type==0||e->m_type==1){
        m_nodespq.push(e);
    }else{
        m_mpq.push(e);
    }
}


void PartsStore::Parts::putSBB(xSBB& b) {
    m_sbbs.emplace_back(b);
    double s=max(m_ps->m_simpquery.m_startTime(), b.m_startTime);
    double e=min(m_ps->m_simpquery.m_endTime(), b.m_endTime);
    auto it = m_line.begin();
    while(it != m_line.end() && (it->ts >= e || it->te <= s)) it++;
    bool isfirst = b.m_startTime<=m_ps->m_simpquery.m_startTime(),
            islast = b.m_endTime>= m_ps->m_simpquery.m_endTime();

    if(isfirst) m_firstsbb = b;
    if(islast) m_lastsbb = b;
    if(it->ts==s&&it->te==e){
        it->d=m_ps->m_simpquery.sbbDist(b);
        auto lit=it,nit =it;
        lit--;nit++;
        if(lit != m_line.end() && !lit->d.infer){
            it->ts=lit->ts;
            it->d+=lit->d;
            m_line.erase(lit);
        }
        if(nit != m_line.end() && !nit->d.infer){
            it->te=nit->te;
            it->d+=nit->d;
            m_line.erase(nit);
        }
        return;
    }
    if(s!=it->ts){ //fill prev one
        double gaps=it->ts,gape=s;
        if(b.m_startTime<=m_mintime&&!m_hasPrev){
            m_line.insert(it, slab(gaps, gape,
                                   DISTE(m_ps->m_simpquery.frontDistStatic(b).opt, 1e300, false)));
            m_firstsbb = b;
        }else{
            m_line.insert(it, slab(gaps, gape,
                                   m_ps->m_simpquery.frontDist(b, tjstat->vmax * 2)));
        }
    }
    if(e!=it->te){ //fill next one
        double gaps=e,gape=it->te;
        auto nit = it;
        nit++;
        if(b.m_endTime>=m_maxtime&&!m_hasNext){
            m_line.insert(nit, slab(gaps, gape,
                                    DISTE(m_ps->m_simpquery.backDistStatic(b).opt, 1e300, false)));
            m_lastsbb = b;
        }else{
            m_line.insert(nit, slab(gaps, gape,
                                    m_ps->m_simpquery.backDist(b, tjstat->vmax * 2)));
        }
    }
    if(s==it->ts){//try to merge into prev one
        auto lit=it;
        lit--;
        if(lit != m_line.end() && !lit->d.infer && lit != m_line.begin()){ //prevent the first
            lit->te= e;
            lit->d+=m_ps->m_simpquery.sbbDist(b);
            m_line.erase(it);
            return;
        }
    }else if(e==it->te){//try to merge into next one
        auto nit=it;
        nit++;
        auto nnit = nit;
        nnit++;
        if(nit != m_line.end() && !nit->d.infer && nnit != m_line.end()){//prevent the last
            nit->ts= s;
            nit->d+=m_ps->m_simpquery.sbbDist(b);
            m_line.erase(it);
            return;
        }
    }
    it->ts=s;
    it->te=e;
    it->d=m_ps->m_simpquery.sbbDist(b);
}