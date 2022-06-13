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
    for(auto &b:parts->m_UCsbbs){
        parts->putSBB(b);
    }
    parts->m_UCsbbs.clear();
    DISTE res;
    int type = 2;
    if (parts->m_missingLeaf.empty()) type = 3;

    if(type == 3) {
        if (parts->m_firstsbb.m_startTime > m_simpquery.m_startTime())
        {
            parts->m_line.front().d.opt = m_simpquery.frontDistStatic(parts->m_firstsbb, parts->m_maxe).opt;
        }
        if (parts->m_lastsbb.m_endTime < m_simpquery.m_endTime())
        {
            parts->m_line.back().d.opt = m_simpquery.backDistStatic(parts->m_lastsbb, parts->m_maxe).opt;
        }
    }

    for (auto it = parts->m_line.begin(); it != parts->m_line.end(); it++) {
        if (!it->d.infer) {
            res += it->d;
        } else {
            if (!m_nodespq.empty())
                res += DISTE(
                        min(it->d.opt, (it->te - it->ts) * m_nodespq.top()->m_dist.opt /
                        (m_simpquery.m_endTime() -
                         m_simpquery.m_startTime())), 1e300, true);
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


void PartsStoreBFMST::insert(id_type id, xTrajectory &r, id_type prev, id_type next, xStoreEntry &entry){
    if(m_parts.find(id)==m_parts.end()){
        m_parts[id]=Parts(this);
    }
//                    m_parts[id].insert(r,prev,next,entry);
    xTrajectory tmp;
    tmp.m_points.emplace_back(r.m_points[0]);
    tmp.m_points.emplace_back(r.m_points[r.m_points.size()-1]);
    m_parts[id].m_computedDist[std::make_pair(r.m_startTime(),r.m_endTime())]
            = DISTE( r.getMinimumDistance(m_query));
    m_parts[id].insert(tmp,prev,next,entry);
    if(global_maxe > m_parts[id].m_maxe)
        m_parts[id].m_maxe = global_maxe;
//                    std::cerr<<"part"<<id<<"\t"<<m_parts[id].m_loadedTime<<"\t"<<m_query.m_endTime()-m_query.m_startTime()<<"\n";
    update(id);
}

DISTE PartsStoreBFMST::update(id_type id) {
    if(m_except.find(id)!=m_except.end()){
        return DISTE(1e300);
    }
    Parts *parts = &m_parts[id];
    double computedTime = 0;
    DISTE pd(0);
    DISTE res;
    std::pair<double, double> timeInterval;


    //mid dist
    const xTrajectory *prev= nullptr;
    for (const auto &traj:parts->m_pTrajs) {
        //this box
        timeInterval.first=traj.m_startTime();
        timeInterval.second=traj.m_endTime();
        if (parts->m_computedDist.count(timeInterval) > 0) {
            if (!parts->m_computedDist[timeInterval].infer) {
                pd = parts->m_computedDist[timeInterval];
            } else {
                pd = DISTE(traj.getMinimumDistance(m_query));
                parts->m_computedDist[timeInterval] = pd;
                if(global_maxe > parts->m_maxe)
                    parts->m_maxe = global_maxe;
            }
        } else {
            pd = DISTE(traj.getMinimumDistance(m_query));
            parts->m_computedDist[timeInterval] = pd;
            if(global_maxe > parts->m_maxe)
                parts->m_maxe = global_maxe;
        }
        res = res+pd;
        computedTime += timeInterval.second - timeInterval.first;
        //the gap
        if (traj.m_startTime() != parts->m_pTrajs.front().m_startTime()) {//not first
            if (prev->m_endTime() < traj.m_startTime()) {
                timeInterval.first=prev->m_endTime();
                timeInterval.second=traj.m_startTime();
                if (parts->m_computedDist.count(timeInterval) > 0) {
                    pd= parts->m_computedDist[timeInterval];
                } else {
                    pd = m_query.gapDist(prev->m_points.back(), traj.m_points.front(), stat->vmax);
                    parts->m_computedDist[timeInterval] = pd;
                }
                if(parts->m_computedDist[timeInterval].infer&&!m_nodespq.empty()) {
                    pd.opt = std::max(pd.opt, (m_lastNodeDist *
                                              (timeInterval.second - timeInterval.first) /
                                              (m_query.m_endTime() - m_query.m_startTime())));
                    pd.pes = max(pd.opt,pd.pes);
                }
                res = res+pd;
            }
        }
        prev = &traj;
    }
    //inferred distance(front dist, back dist and mid dist) should be stored as negative values
    //front dist
    if (parts->m_mintime > m_query.m_startTime()) {
        timeInterval.first=m_query.m_startTime();
        timeInterval.second=parts->m_mintime;
        if (parts->m_computedDist.count(timeInterval) > 0) {
            pd= parts->m_computedDist[timeInterval];
        } else {
            if (parts->m_hasPrev) {
                pd = m_query.frontDist(parts->m_pTrajs.front().m_points.front(),
                                       stat->vmax);
                parts->m_computedDist[timeInterval] = pd;
            } else {
                pd = m_query.frontDistStatic(parts->m_pTrajs.front().m_points.front(),parts->m_maxe);
                parts->m_computedDist[timeInterval] = pd;
                computedTime += timeInterval.second - timeInterval.first;
            }
        }
        if (parts->m_computedDist[timeInterval].infer &&
            !m_nodespq.empty()) {
            pd.opt = std::max(pd.opt, m_lastNodeDist *
                                      (timeInterval.second -
                                       timeInterval.first) /
                                      (m_query.m_endTime() -
                                       m_query.m_startTime()));
            pd.pes = max(pd.opt, pd.pes);
        }
        res = res+pd;
    }
    //backdist
    if (parts->m_maxtime < m_query.m_endTime()) {
        timeInterval.first=parts->m_maxtime;
        timeInterval.second=m_query.m_endTime();
        if (parts->m_computedDist.count(timeInterval) > 0) {
            pd= parts->m_computedDist[timeInterval];
        } else {
            if (parts->m_hasNext) {
                pd = m_query.backDist(parts->m_pTrajs.back().m_points.back(),stat->vmax);
                parts->m_computedDist[timeInterval] = pd;
            } else {
                pd = m_query.backDistStatic(parts->m_pTrajs.back().m_points.back(),parts->m_maxe);
                parts->m_computedDist[timeInterval] = pd;
                computedTime += timeInterval.second - timeInterval.first;
            }
        }
            if (parts->m_computedDist[timeInterval].infer &&
                !m_nodespq.empty()) {
                pd.opt = std::max(pd.opt, m_lastNodeDist *
                                          (timeInterval.second -
                                           timeInterval.first) /
                                          (m_query.m_endTime() -
                                           m_query.m_startTime()));
                pd.pes = max(pd.opt, pd.pes);
            }
        res = res+pd;
    }

    parts->m_calcMin = res;
    parts->m_computedTime = computedTime;
    int type = 2;
    if (parts->m_loadedTime + 1e-7 >= (m_query.m_endTime() - m_query.m_startTime())) type = 3;
    res.opt -= m_error;
    res.pes += m_error;
    if (m_handlers.find(id) == m_handlers.end()) {
        auto handle = m_mpq.push(new NNEntry(id, res, type));
        m_handlers[id] = handle;
    } else {
        m_mpq.update(m_handlers[id], id, res, type);
    }
    m_pes.insert(id, res.pes);
    if(res.opt>m_pes.threshold()){
        m_except.insert(id);
    }
    return res;
}

void PartsStoreBFMST::loadPartTraj(id_type id, leafInfo *e, double dist){
    xTrajectory tmpTraj;
    id_type trajid = e->m_se.m_id;
    m_ts->loadTraj(tmpTraj, e->m_se);
    if (m_except.find(trajid) != m_except.end()) return;
    xTrajectory inter;
    tmpTraj.getPartialxTrajectory(max(m_query.m_startTime(), e->b.m_startTime),
                                  min(m_query.m_endTime(), e->b.m_endTime), inter);
#ifdef TJDEBUG
//    cerr<<"load traj into id "<<trajid<<endl<<inter.toString()<<endl;
#endif
    bool pv = e->m_hasPrev, nt = e->m_hasNext;
    if (!m_pTree->m_bStoringLinks) {
        pv = e->m_se.m_s > 0 || e->b.m_startTime > tmpTraj.m_startTime();
        xTrajEntry xte = db_load_traj_entry(e->m_se.m_id);
        nt = e->m_se.m_e < xte.m_npoint-1 || e->b.m_endTime < tmpTraj.m_endTime();
    }
    pv = pv && (m_query.m_startTime() < e->b.m_startTime);
    nt = nt && (m_query.m_endTime() > e->b.m_endTime);
    insert(trajid, inter, pv ? 1 : -1, nt ? 1 : -1, e->m_se);
}

NNEntry* PartsStoreBFMST::top() {
    if(m_mpq.empty()||(!m_nodespq.empty()&&m_nodespq.top()->m_dist< m_mpq.top()->m_dist)){
        return m_nodespq.top();
    }
    if(!m_mpq.empty()){
        id_type prevtop=-1;
        while(m_mpq.top()->m_id!=prevtop){
            while(m_except.find(m_mpq.top()->m_id)!=m_except.end()){
                auto s= m_mpq.top();
                m_mpq.pop();
                delete s;
                if(m_mpq.empty()){
                    return m_nodespq.top();
                }
            }
            prevtop=m_mpq.top()->m_id;
            update(prevtop);
        }
    }
    if(!m_mpq.empty()&&m_mpq.top()->m_type==3&&(m_nodespq.empty()|| m_mpq.top()->m_dist < m_nodespq.top()->m_dist)){
        return m_mpq.top();
    }
    if(!m_nodespq.empty())
        return m_nodespq.top();
    return m_mpq.top();
}

NNEntry* PartsStoreBFMST::pop(){
    if(!m_mpq.empty()&&m_mpq.top()->m_type==3&&(m_nodespq.empty()|| m_mpq.top()->m_dist < m_nodespq.top()->m_dist))
        return m_mpq.pop();
    if(!m_nodespq.empty()) {
        if(m_nodespq.top()->m_type==0){
            m_lastNodeDist = m_nodespq.top()->m_dist.opt;
        }
        return m_nodespq.pop();
    }
    return m_mpq.pop();
}

string PartsStoreBFMST::explain(id_type id){
    if(m_parts.find(id)==m_parts.end()) return "";
     stringstream  ss;
     auto s =m_parts[id];
     ss<<"id is "<<id<<endl;
     for(auto &b:m_parts[id].m_pTrajs){
         ss<<b.toString()<<endl;
     }
     for(auto &b:m_parts[id].m_computedDist){
         if(b.second.infer== false){
             ss<<b.first.first<<"\t"<<b.first.second<<"\t"<<b.second.opt<<endl;
         }
     }
     return ss.str();
 }

void PartsStoreBFMST::push(NNEntry *e) {
    if(e->m_type==0||e->m_type==1){
        m_nodespq.push(e);
    }else{
        m_mpq.push(e);
    }
}

void PartsStore::Parts::putSBB(xSBB& b) {
    #ifndef NDEBUG
    m_sbbs.emplace_back(b);
    #endif
    double s=max(m_ps->m_simpquery.m_startTime(), b.m_startTime);
    double e=min(m_ps->m_simpquery.m_endTime(), b.m_endTime);
    auto it = m_line.begin();
    while(it != m_line.end() && (it->ts >= e || it->te <= s)) it++;
    if(b.m_startTime<m_ps->m_simpquery.m_startTime()) m_firstsbb = b;
    if(b.m_endTime> m_ps->m_simpquery.m_endTime()) m_lastsbb = b;
    if(it->ts==s&&it->te==e){
        it->d=m_ps->m_simpquery.sbbDist(b);
        if(global_maxe>m_maxe)
            m_maxe=global_maxe;
        auto lit=it,nit =it;
        lit--;nit++;
        if(lit != m_line.end() && !lit->d.infer){
            it->ts=lit->ts;
            it->d+=lit->d;
            m_line.erase(lit);
        }
        if(nit != m_line.end() && !nit->d.infer){
            it->ts=nit->ts;
            it->d+=nit->d;
            m_line.erase(nit);
        }
        return;
    }
    if(s!=it->ts){ //fill prev one
        double gaps=it->ts,gape=s;
        if(b.m_startTime<=m_mintime&&!m_hasPrev){
            m_line.insert(it, slab(gaps, gape,
                                   DISTE(m_ps->m_simpquery.frontDistStatic(b, it->d.opt/(gape-gaps)).opt,1e300, false)));
            m_firstsbb = b;
        }else{
            m_line.insert(it, slab(gaps, gape,
                                   m_ps->m_simpquery.frontDist(b, tjstat->vmax*2)));
        }
    }
    if(e!=it->te){ //fill next one
        double gaps=e,gape=it->te;
        auto nit = it;
        nit++;
        if(b.m_endTime>=m_maxtime&&!m_hasNext){
            m_line.insert(nit, slab(gaps, gape,
                                    DISTE(m_ps->m_simpquery.backDistStatic(b, it->d.opt/(gape-gaps)).opt,1e300, false)));
            m_lastsbb = b;
        }else{
            m_line.insert(nit, slab(gaps, gape,
                                   m_ps->m_simpquery.backDist(b, tjstat->vmax*2)));
        }
    }
    if(s==it->ts){//try to merge into prev one
        auto lit=it;
        lit--;
        if(lit != m_line.end() && !lit->d.infer && lit != m_line.begin()){ //prevent the first
            lit->te= e;
            lit->d+=m_ps->m_simpquery.sbbDist(b);
            if(global_maxe>m_maxe)
                m_maxe=global_maxe;
            m_line.erase(it);
            return;
        }
    }else if(e==it->te){//try to merge into next one
        auto nit=it;
        nit++;
        auto nnit = nit++;
        if(nit != m_line.end() && !nit->d.infer && nnit != m_line.end()){//prevent the last
            nit->ts= s;
            nit->d+=m_ps->m_simpquery.sbbDist(b);
            if(global_maxe>m_maxe)
                m_maxe=global_maxe;
            m_line.erase(it);
            return;
        }
    }
    it->ts=s;
    it->te=e;
    it->d=m_ps->m_simpquery.sbbDist(b);
    if(global_maxe>m_maxe)
        m_maxe=global_maxe;
}