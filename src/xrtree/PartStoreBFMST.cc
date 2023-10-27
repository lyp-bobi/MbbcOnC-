//
// Created by Chuang on 2023/10/23.
//

#include "xRTree.h"

using namespace xRTreeNsp;

void PartsStoreBFMST::insert(id_type id, xTrajectory &r, id_type prev, id_type next, xStoreEntry &entry) {
    if (m_parts.find(id) == m_parts.end()) {
        m_parts[id] = Parts(this);
    }
//                    m_parts[id].insert(r,prev,next,entry);
    if (current_distance != IED) {
        m_parts[id].m_computedDist[std::make_pair(r.m_startTime(),
                                                  r.m_endTime())]
                = m_query.getPartialRMDTW(r);
        m_parts[id].insert(r, prev, next, entry);
    } else{
        xTrajectory tmp;
        tmp.m_points.emplace_back(r.m_points[0]);
        tmp.m_points.emplace_back(r.m_points[r.m_points.size()-1]);
        m_parts[id].m_computedDist[std::make_pair(r.m_startTime(),r.m_endTime())]
                = DISTE( r.getMinimumDistance(m_query));
        m_parts[id].insert(tmp,prev,next,entry);
    }
//                    std::cerr<<"part"<<id<<"\t"<<m_parts[id].m_loadedTime<<"\t"<<m_query.m_endTime()-m_query.m_startTime()<<"\n";
    update(id);
}

DISTE PartsStoreBFMST::update(id_type id) {
    Parts *parts = &m_parts[id];
    double computedTime = 0;
    DISTE pd(0);
    DISTE res;
    std::pair<double, double> timeInterval;

    int type = 2;
    if (parts->m_loadedTime + 1e-7 >= (m_query.m_endTime() - m_query.m_startTime())) type = 3;

    if(type == 3 && current_distance != IED)
    {
        xTrajectory fulltraj;
        for(auto &pj: parts->m_pTrajs)
        {
            if(fulltraj.m_points.size() == 0)
                fulltraj = pj;
            else
                fulltraj.linkxTrajectory(pj);
        }
        res = DISTE(m_query.getMinimumDistance(fulltraj));
        parts->m_calcMin = res;
        parts->m_computedTime = parts->m_loadedTime;
    } else {


        //mid dist
        const xTrajectory *prev = nullptr;
        for (const auto &traj:parts->m_pTrajs) {
            //this box
            timeInterval.first = traj.m_startTime();
            timeInterval.second = traj.m_endTime();
            if (parts->m_computedDist.count(timeInterval) > 0) {
                if (!parts->m_computedDist[timeInterval].infer) {
                    pd = parts->m_computedDist[timeInterval];
                } else {
                    if (current_distance != IED)
                        pd = m_query.getPartialRMDTW(traj);
                    else
                        pd = DISTE(traj.getMinimumDistance(m_query));
                    parts->m_computedDist[timeInterval] = pd;
                }
            } else {
                if (current_distance != IED)
                    pd = m_query.getPartialRMDTW(traj);
                else
                    pd = DISTE(traj.getMinimumDistance(m_query));
                parts->m_computedDist[timeInterval] = pd;
            }
            res = res + pd;
            computedTime += timeInterval.second - timeInterval.first;
            //the gap
            if (current_distance == IED) {
                if (traj.m_startTime() !=
                    parts->m_pTrajs.front().m_startTime()) {//not first
                    if (prev->m_endTime() < traj.m_startTime()) {
                        timeInterval.first = prev->m_endTime();
                        timeInterval.second = traj.m_startTime();
                        if (parts->m_computedDist.count(timeInterval) > 0) {
                            pd = parts->m_computedDist[timeInterval];
                        } else {
                            pd = m_query.gapDist(prev->m_points.back(),
                                                 traj.m_points.front(),
                                                 stat->vmax);
                            parts->m_computedDist[timeInterval] = pd;
                        }
                        if (parts->m_computedDist[timeInterval].infer &&
                            !m_nodespq.empty()) {
                            pd.opt = std::max(pd.opt, (m_lastNodeDist *
                                                       (timeInterval.second -
                                                        timeInterval.first) /
                                                       (m_query.m_endTime() -
                                                        m_query.m_startTime())));
                            pd.pes = max(pd.opt, pd.pes);
                        }
                        res = res + pd;
                    }
                }
                prev = &traj;
            }
        }
        //inferred distance(front dist, back dist and mid dist) should be stored as negative values
        //front dist
        if (parts->m_mintime > m_query.m_startTime()) {
            pd = DISTE(0);
            timeInterval.first = m_query.m_startTime();
            timeInterval.second = parts->m_mintime;
            if (parts->m_computedDist.count(timeInterval) > 0) {
                pd = parts->m_computedDist[timeInterval];
            } else {
                if (parts->m_hasPrev) {
//                pd = m_query.frontDist(parts->m_pTrajs.front().m_points.front(),
//                                       stat->vmax);
//                parts->m_computedDist[timeInterval] = pd;
                } else {
                    pd = m_query.frontDistStatic(
                            parts->m_pTrajs.front().m_points.front());
                    parts->m_computedDist[timeInterval] = pd;
                    computedTime += timeInterval.second - timeInterval.first;
                }
            }
//        if (parts->m_computedDist[timeInterval].infer &&
//            !m_nodespq.empty()) {
//            pd.opt = std::max(pd.opt, m_lastNodeDist *
//                                      (timeInterval.second -
//                                       timeInterval.first) /
//                                      (m_query.m_endTime() -
//                                       m_query.m_startTime()));
//            pd.pes = max(pd.opt, pd.pes);
//        }
            res = res + pd;
        }
        //backdist
        if (parts->m_maxtime < m_query.m_endTime()) {
            pd = DISTE(0);
            timeInterval.first = parts->m_maxtime;
            timeInterval.second = m_query.m_endTime();
            if (parts->m_computedDist.count(timeInterval) > 0) {
                pd = parts->m_computedDist[timeInterval];
            } else {
                if (parts->m_hasNext) {
//                pd = m_query.backDist(parts->m_pTrajs.back().m_points.back(),stat->vmax);
//                parts->m_computedDist[timeInterval] = pd;
                } else {
                    pd = m_query.backDistStatic(
                            parts->m_pTrajs.back().m_points.back());
                    parts->m_computedDist[timeInterval] = pd;
                    computedTime += timeInterval.second - timeInterval.first;
                }
            }
//            if (parts->m_computedDist[timeInterval].infer &&
//                !m_nodespq.empty()) {
//                pd.opt = std::max(pd.opt, m_lastNodeDist *
//                                          (timeInterval.second -
//                                           timeInterval.first) /
//                                          (m_query.m_endTime() -
//                                           m_query.m_startTime()));
//                pd.pes = max(pd.opt, pd.pes);
//            }
            res = res + pd;
        }
    }
    res.opt -= m_error;
    res.pes += m_error;
    parts->m_calcMin = res;
    parts->m_computedTime = computedTime;
    if (m_handlers.find(id) == m_handlers.end()) {
        auto handle = m_mpq.push(new NNEntry(id, res, type));
        m_handlers[id] = handle;
    } else {
        m_mpq.update(m_handlers[id], id, res, type);
    }
    return res;
}

void PartsStoreBFMST::loadPartTraj(id_type id, leafInfo *e, double dist){
    xTrajectory tmpTraj;
    id_type trajid = e->m_se.m_id;
    m_ts->loadTraj(tmpTraj, e->m_se);
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
//         ss<<b.toString()<<endl;
        ss << b.m_startTime() << "\t"<< b.m_endTime()<<endl;
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
