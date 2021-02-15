//
// Created by Chuang on 2021/2/3.
//

#include "xRTree.h"

using namespace xRTreeNsp;

void PartsStore::Parts::insert(xSBB &r, id_type prev, id_type next, xStoreEntry &entry) {
    is_modified = true;
    m_sbbs[r.m_startTime]=r;
    m_loadedTime += min(r.m_endTime, m_ps->m_query.m_endTime()) -
                    max(r.m_startTime, m_ps->m_query.m_startTime());
    if(r.m_endTime > m_maxtime){
        m_maxtime=min(r.m_endTime, m_ps->m_query.m_endTime());
        if(next==-1|| r.m_endTime > m_ps->m_query.m_endTime()) {
            m_hasNext=false;
            m_loadedTime+=m_ps->m_query.m_endTime()-m_maxtime;
        }
    }
    if(r.m_startTime < m_mintime){
        m_mintime=max(r.m_startTime, m_ps->m_query.m_startTime());
        if(prev==-1|| r.m_startTime < m_ps->m_query.m_startTime()) {
            m_hasPrev=false;
            m_loadedTime+=m_mintime-m_ps->m_query.m_startTime();
        }
    }
    if(prev>=0&&m_loadedLeaf.count(prev)==0) m_missingLeaf.insert(prev);
    if(next>=0&&m_loadedLeaf.count(next)==0) m_missingLeaf.insert(next);
    m_ses[r.m_startTime]=entry;
}

void PartsStore::insert(id_type id, xSBB &b, id_type prev, id_type next, xStoreEntry &entry) {
    if(m_parts.count(id)==0){
        m_parts[id]=Parts(this);
    }
    m_parts[id].insert(b,prev,next,entry);
}

DISTE PartsStore::updateValue(id_type id) {
    Parts *parts = &m_parts[id];
    if(!parts->is_modified) return parts->m_calcMin;
    double computedTime = 0;
    DISTE pd;
    DISTE res;
    bool needComp =false;
    std::pair<double, double> timeInterval;
    std::map<std::pair<double,double>,DISTE>::iterator  iter;
    //inferred distance(front dist, back dist and mid dist) should be stored as negative values
    //front dist
    if (parts->m_mintime > m_query.m_startTime()) {
        timeInterval.first=m_query.m_startTime();
        timeInterval.second=parts->m_mintime;
        iter = parts->m_computedDist.find(timeInterval);
        if(iter != parts->m_computedDist.end()){
            pd= iter->second;
        } else {
            if (parts->m_hasPrev) {
                pd = m_query.frontDist(parts->m_sbbs.begin()->second,stat->vmax);
                parts->m_computedDist[timeInterval] = pd;
            } else {
                pd = m_query.frontDistStatic(parts->m_sbbs.begin()->second);
                parts->m_computedDist[timeInterval] = pd;
                computedTime += timeInterval.second - timeInterval.first;
            }
        }
//        //we don't know if it has sbbs, so can't use
//        if(tjstat->regular&&iter->second.infer&&!m_nodespq.empty()){
//            pd.opt = std::max(pd.opt, m_nodespq.top()->m_dist.opt *
//                                      (timeInterval.second - timeInterval.first) /
//                                      (m_query.m_endTime() - m_query.m_startTime()));
//            pd.pes = max(pd.opt,pd.pes);
//        }
        res = res + pd;
    }
    //mid dist
    const xSBB *prev= nullptr;
    for (const auto &box:parts->m_sbbs) {
        //this box
        timeInterval.first=max(box.second.m_startTime, m_query.m_startTime());
        timeInterval.second=min(box.second.m_endTime, m_query.m_endTime());
        iter = parts->m_computedDist.find(timeInterval);
        if(iter != parts->m_computedDist.end()&&!(iter->second.infer)){
            pd= iter->second;
        } else {
            pd = DISTE(m_query.sbbDist(box.second));
            parts->m_computedDist[timeInterval] = pd;
        }
        res = res + pd;
        computedTime += timeInterval.second - timeInterval.first;
        //the gap
        if (box.second.m_startTime != parts->m_sbbs.begin()->second.m_startTime) {//not first
            if (prev->m_endTime < box.second.m_startTime) {
                timeInterval.first= prev->m_endTime;
                timeInterval.second= box.second.m_startTime;
                iter = parts->m_computedDist.find(timeInterval);
                if(iter != parts->m_computedDist.end()){
                    pd= iter->second;
                } else {
                    pd = m_query.gapDist(*prev, box.second, stat->vmax);
                    parts->m_computedDist[timeInterval] = pd;
                }
                if(iter->second.infer&&!m_nodespq.empty()){
                    pd.opt = std::max(pd.opt, ( m_nodespq.top()->m_dist.opt *
                                              (timeInterval.second - timeInterval.first) /
                                              (m_query.m_endTime() - m_query.m_startTime())));
                    pd.pes = max(pd.opt,pd.pes);
                }
                res = res + pd;
            }
        }
        prev = &box.second;
    }
    //backdist
    if (parts->m_maxtime < m_query.m_endTime()) {
        timeInterval.first=parts->m_maxtime;
        timeInterval.second=m_query.m_endTime();
        iter = parts->m_computedDist.find(timeInterval);
        if(iter != parts->m_computedDist.end()){
            pd= iter->second;
        } else {
            if (parts->m_hasNext) {
                pd = m_query.backDist(parts->m_sbbs.rbegin()->second,stat->vmax);
                parts->m_computedDist[timeInterval] = pd;
            } else {
                pd = m_query.backDistStatic(parts->m_sbbs.rbegin()->second);
                parts->m_computedDist[timeInterval] = pd;
                computedTime += timeInterval.second - timeInterval.first;
            }
        }
//        //we don't know if it has sbbs, so can't use
//        if(tjstat->regular&&iter->second.infer&&!m_nodespq.empty()){
//            pd.opt = std::max(pd.opt, m_nodespq.top()->m_dist.opt *
//                                      (timeInterval.second - timeInterval.first) /
//                                      (m_query.m_endTime() - m_query.m_startTime()));
//            pd.pes = max(pd.opt,pd.pes);
//        }
        res = res + pd;
    }
    parts->m_calcMin = res;
    parts->m_computedTime = computedTime;
    int type = 2;
    if (parts->m_missingLeaf.empty()) type = 3;
    res.opt -= m_error;
    res.pes += m_error;
    m_pes.insert(id, res.pes);
    if(res.opt>m_pes.threshold()){
        m_except.insert(id);
    }
    m_mpq.updateValue(m_handlers[id], id, res, type);
    parts->is_modified=false;
    return res;
}


void PartsStore::loadLeaf(const Node &n, double dist) {
//                    std::cerr<<"load leaf"<<n.m_nodeMBR<<"\n";
//                    std::cerr<<"leaf dist"<<m_query.getNodeMinimumDistance(n.m_nodeMBR,100)/(m_query.m_endTime()-m_query.m_startTime())<<"\n";
//                    std::cerr<<"load leaf"<<n.m_identifier<<"\n";
    loadedLeaf.insert(n.m_identifier);
    std::set<id_type > relatedIds;
    for(int i=0;i<n.m_children;i++){
        id_type trajid=n.m_se[i].m_id;
        xStoreEntry entry= n.m_se[i];
        double bts= n.m_ptrxSBB[i]->m_startTime,bte= n.m_ptrxSBB[i]->m_endTime;
        if(bts>=m_query.m_endTime()||
           bte<=m_query.m_startTime()||m_except.count(trajid)>0){}
        else {
            insert(trajid, *n.m_ptrxSBB[i],
                   (m_query.m_startTime()<bts)?n.m_prevNode[i]:-1
                    , (m_query.m_endTime()>bte)?n.m_nextNode[i]:-1,
                   entry);
            relatedIds.insert(trajid);
//                        if(m_parts.count(trajid)>0){
//                            auto s = m_parts[trajid];
//                            std::cerr<<"";
//                        }
        }
    }
    for(const auto &rid:relatedIds){
        m_parts[rid].m_missingLeaf.erase(n.m_identifier);
        m_parts[rid].m_loadedLeaf.insert(n.m_identifier);
        if(m_handlers.count(rid)==0){
            auto handle = m_mpq.push(new NNEntry(rid, DISTE(dist), 2));
            m_handlers[rid] = handle;
        }
    }
}

NNEntry* PartsStore::top() {
    id_type lastid = -1;
    if(m_nodespq.empty()) {
        while(m_mpq.top()->m_type==2&&lastid!=m_mpq.top()->m_id) {
            lastid = m_mpq.top()->m_id;
            updateValue(lastid);
            m_mpq.updateOrder(m_handlers[lastid]);
        }
        return m_mpq.top();
    }
    while(true){
        if(m_mpq.empty()||m_nodespq.top()->m_dist < m_mpq.top()->m_dist){
            return m_nodespq.top();
        }
        if(m_mpq.top()->m_type==4)
            return m_mpq.top();
        if(lastid!=m_mpq.top()->m_id) {
            lastid = m_mpq.top()->m_id;
            updateValue(lastid);
            m_mpq.updateOrder(m_handlers[lastid]);
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
    if(m_parts.count(id)==0){
        m_parts[id]=Parts(this);
    }
//                    m_parts[id].insert(r,prev,next,entry);
    xTrajectory tmp;
    tmp.m_points.emplace_back(r.m_points[0]);
    tmp.m_points.emplace_back(r.m_points[r.m_points.size()-1]);
    m_parts[id].m_computedDist[std::make_pair(r.m_startTime(),r.m_endTime())]
            = DISTE( r.getMinimumDistance(m_query));
    m_parts[id].insert(tmp,prev,next,entry);
//                    std::cerr<<"part"<<id<<"\t"<<m_parts[id].m_loadedTime<<"\t"<<m_query.m_endTime()-m_query.m_startTime()<<"\n";
    update(id);
}

DISTE PartsStoreBFMST::update(id_type id) {
    if(m_except.count(id)>0){
        return DISTE(1e300);
    }
    Parts *parts = &m_parts[id];
    double computedTime = 0;
    DISTE pd(0);
    DISTE res;
    std::pair<double, double> timeInterval;

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
                pd = m_query.frontDistStatic(parts->m_pTrajs.front().m_points.front());
                parts->m_computedDist[timeInterval] = pd;
                computedTime += timeInterval.second - timeInterval.first;
            }
        }
        //we don't know if it has sbbs, so can't use
//                    if(parts->m_computedDist[timeInterval].infer&&!m_nodespq.empty()){
//                        pd.opt = std::max(pd.opt, m_lastNodeDist *
//                                                  (timeInterval.second - timeInterval.first) /
//                                                  (m_query.m_endTime() - m_query.m_startTime()));
//                        pd.pes = max(pd.opt,pd.pes);
//                    }
        res = res+pd;
    }
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
            }
        } else {
            pd = DISTE(traj.getMinimumDistance(m_query));
            parts->m_computedDist[timeInterval] = pd;
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
                pd = m_query.backDistStatic(parts->m_pTrajs.back().m_points.back());
                parts->m_computedDist[timeInterval] = pd;
                computedTime += timeInterval.second - timeInterval.first;
            }
        }
        //we don't know if it has sbbs, so can't use
//                    if(parts->m_computedDist[timeInterval].infer&&!m_nodespq.empty()){
//                        pd.opt = std::max(pd.opt, m_lastNodeDist *
//                                                  (timeInterval.second - timeInterval.first) /
//                                                  (m_query.m_endTime() - m_query.m_startTime()));
//                        pd.pes = max(pd.opt,pd.pes);
//                    }
        res = res+pd;
    }

    parts->m_calcMin = res;
    parts->m_computedTime = computedTime;
    int type = 2;
    if (parts->m_loadedTime + 1e-7 >= (m_query.m_endTime() - m_query.m_startTime())) type = 3;
    res.opt -= m_error;
    res.pes += m_error;
    if (m_handlers.count(id) == 0) {
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
    if (m_except.count(trajid) > 0) return;
    xTrajectory inter;
    tmpTraj.getPartialxTrajectory(max(m_query.m_startTime(), e->m_ts),
                                  min(m_query.m_endTime(), e->m_te), inter);
#ifdef TJDEBUG
//    cerr<<"load traj into id "<<trajid<<endl<<inter.toString()<<endl;
#endif
    bool pv = e->m_hasPrev, nt = e->m_hasNext;
    if (!m_pTree->m_bStoringLinks) {
        pv = e->m_se.m_s > 0 || e->m_ts > tmpTraj.m_startTime();
        nt = e->m_se.m_e < m_ts->m_trajIdx[e->m_se.m_id]->m_npoint-1 || e->m_te < tmpTraj.m_endTime();
    }
    pv = pv && (m_query.m_startTime() < e->m_ts);
    nt = nt && (m_query.m_endTime() > e->m_te);
    insert(trajid, inter, pv ? 1 : -1, nt ? 1 : -1, e->m_se);
}

NNEntry* PartsStoreBFMST::top() {
    if(m_mpq.empty()||(!m_nodespq.empty()&&m_nodespq.top()->m_dist< m_mpq.top()->m_dist)){
        return m_nodespq.top();
    }
    if(!m_mpq.empty()){
        id_type prevtop=-1;
        while(m_mpq.top()->m_id!=prevtop){
            while(m_except.count(m_mpq.top()->m_id)>0){
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
    if(m_parts.count(id)==0) return "";
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