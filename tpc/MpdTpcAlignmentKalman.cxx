#include "MpdTpcAlignmentKalman.h"
#include <TVector3.h>

MpdTpcAlignmentKalman::StatBuffer& MpdTpcAlignmentKalman::StatBuffer::init(size_t size)
{
    erase();
    size_max = size;
    data = new Eigen::Matrix<double, npars, 1>[size_max];

    return *this;
}


MpdTpcAlignmentKalman::StatBuffer& MpdTpcAlignmentKalman::StatBuffer::erase()
{
    if (data) delete[] data;
    data = nullptr;

    size_max = offset = 0;
    data_ready = false;

    return *this;
}


MpdTpcAlignmentKalman::StatBuffer& MpdTpcAlignmentKalman::StatBuffer::add(const Eigen::Matrix<double, npars, 1>& val)
{
    if (data && size_max)
    {
        data[offset] = val;
        ++offset;

        if (offset == size_max)
        {
            offset = 0;
            data_ready  = true;
        }
    }

    return *this;
}


MpdTpcAlignmentKalman::StatBuffer& MpdTpcAlignmentKalman::StatBuffer::stats_compute(Eigen::Matrix<double, npars, 1> &mean, Eigen::Matrix<double, npars, 1> &disp)
{
    mean.setZero();
    disp.setZero();

    auto length = (data_ready ? size_max : offset);

    if (length)
    {
        for (size_t i = 0; i < length; ++i)
            mean += data[i];

        auto norm = 1./length;
        mean *= norm;

        if (length > 1)
        {
            for (size_t i = 0; i < length; ++i)
                disp += ((data[i] - mean).array()*(data[i] - mean).array()).matrix();

            norm = 1./(length-1);
            disp *= norm;

        }
    }

    return *this;
}


MpdTpcAlignmentKalman::MpdTpcAlignmentKalman(Option_t* calib_tracks_type)
    : FairTask("TPC Alignment Kalman"), fPersistence(kFALSE)
{
    tracks_type = new char[255];
    std::strcpy(tracks_type, calib_tracks_type);
#ifdef DEBUG
    //reader = nullptr;
    rand_gen = nullptr;
#endif
}


MpdTpcAlignmentKalman::~MpdTpcAlignmentKalman()
{
    delete tracks_type;
}


InitStatus MpdTpcAlignmentKalman::Init()
{
    FairRootManager* frman = FairRootManager::Instance();

    if (!frman)
    {
        Error("MpdTpcAlignmentKalman::Init","RootManager not found");
        return kFATAL;
    }

    if (!(fHits = (TClonesArray*) frman->GetObject("TpcRecPoint")))
    {
        Error("MpdTpcAlignmentKalman::Init","TPC Reco Hits not found");
        return kFATAL;
    }

    hits.reserve(30000);
    tracks.reserve(300);

    // init/fixed values
    Ind.setZero().diagonal().setOnes();
    Inp.setZero().diagonal().setOnes();
    Ins.setZero().diagonal().setOnes();

    A = Ins;

    Q.setZero().diagonal().setConstant(Qini);
    //Q.block(ndims,ndims,ndims,ndims).diagonal().setConstant(0.);

    for (unsigned sid = 0; sid < nSectors; ++sid)
    {
        s[sid].setZero();
        C[sid].setZero().diagonal().setConstant(Cini);

        matrixHInit(sid);
        chi2[sid] = 0.;

        track_segments[sid] = 0;
        ann[sid] = ann0;

        is_finished[sid] = false;
        is_active[sid] = true;

        pstats[sid].init(Awnd);
    }

    exec_time = 0;
    events_processed = tracks_processed = 0;

#ifdef DEBUG
    if (!(fSecGeo = MpdTpcSectorGeo::Instance()))
    {
        Error("MpdTpcMisalignmentKalman::Init","TPC Sector Geo init error");
        return kFATAL;
    }

    /* reader = new MpdTpcAlignmentParams();
    if (reader->ProcessParamsFile())
    {
        Error("MpdTpcAlignmentKalman::Init","Debug params read error");
        return kFATAL;
    } */

    fHitsAlignTrue = new TClonesArray("MpdTpcHit");
    frman->Register("TpcRecPointAlignTrue", "TpcAlign", fHitsAlignTrue, fPersistence);

    out_file = TFile::Open("output.root","RECREATE");
    out_tree = new TTree("tree", "Kalman alignment sums");
    out_tree->Branch("st", ost.t, "t[24]/D:ann[24]/D:nsg[24]/D:nht[24]/D:dx[24]/D:dy[24]/D:dz[24]/D:dp[24]/D:de[24]/D:dk[24]/D:Cx[24]/D:Cy[24]/D:Cz[24]/D:Cp[24]/D:Ce[24]/D:Ck[24]/D");//:Vx[24]/D:Vy[24]/D:Vz[24]/D");

    hits_per_trk = new TH1I("hits_per_trk", "hits_per_trk", 170, 1, 170);
    hits_per_seg = new TH1I("hits_per_seg", "hits_per_seg", 70, 1, 70);
    segs_per_trk = new TH1I("segs_per_trk", "segs_per_trk", 10, 1, 10);

    std::random_device device;
    rand_gen = new std::mt19937(device());
    char buf[255];

    auto p_min = -.7, p_max = .7;

    for (unsigned sid = 0; sid < nSectors; ++sid)
    {
        //st[sid] << 0., 0., reader->GetSectorsShift().fTpcSectorShift[sid].X(), reader->GetSectorsShift().fTpcSectorShift[sid].Y();

        for (int p = 0; p < npars; ++p)
        {
            st[sid](p,0) = 0.;

//            if (!(sid%2))
//                st[sid](p,0) = rnd_flat(p_min, p_max);
//            else
//                is_active[sid] = false;
//            if (!(sid%4))
//                is_active[sid] = false;
//            else
                st[sid](p,0) = rnd_flat(p_min, p_max);
        }

        sprintf(buf, "segs_per_sec_%02d", sid);
        segs_per_sec[sid] = new TH1I(buf, buf, 80, 1, 80);

        sprintf(buf, "hits_per_sec_%02d", sid);
        hits_per_sec[sid] = new TH1I(buf, buf, 1800, 1, 1800);

        //std::printf("Sector %02d: % .7f % .7f % .7f\n", sid, dt[sid](0,0), dt[sid](1,0), dt[sid](2,0));
    }

    weak_modes_fix(st);

    printout_period = 100;
    event_log_flush();
#endif

    Info("MpdTpcAlignmentKalman::Init", "finished OK");
    return kSUCCESS;
}


void MpdTpcAlignmentKalman::Exec(Option_t* opt)
{
    {
        bool is_over =  true;

        for (unsigned sid = 0; sid < nSectors; ++sid)
            is_over = (is_finished[sid] && is_over);

        if (is_over)
            return;
    }

    auto tStart = clock();
    auto nHits  = unsigned(fHits->GetEntriesFast());

    if (!nHits)
    {
        Warning("MpdTpcAlignmentKalman::Exec", "finished with ZERO hits");
        return;
    }

    hits.clear();
    hits.reserve(nHits);

    MpdTpcHit *tpcHit;

#ifdef DEBUG
    misalignment_commit(nHits);

    fHitsAlignTrue->Delete();
    TVector3 locPos, globPos;
#endif

    for (unsigned i = 0; i < nHits; ++i)
    {
        // collecting hits
        tpcHit = static_cast<MpdTpcHit*>(fHits->UncheckedAt(i));
        hits.push_back(hit_(tpcHit));

#ifdef DEBUG
        auto glob = align(hits.back()->sid, hits.back()->loc, st);

        globPos.SetXYZ(glob(0,0), glob(1,0), glob(2,0));
        fSecGeo->Global2Local(globPos, locPos);

        Int_t iHitAlignTrue = fHitsAlignTrue->GetEntriesFast();
        MpdTpcHit *hitAlignTrue = new ((*fHitsAlignTrue)[iHitAlignTrue]) MpdTpcHit(*tpcHit);

        hitAlignTrue->SetPosition(globPos);
        hitAlignTrue->SetLocalPosition(locPos);
#endif
    }

    // finding tracks
    tracks.clear();

    if (tracksFind())
    {
#ifdef DEBUG
        if ( !(events_processed%printout_period))
        {
            std::cout << "MpdTpcAlignmentKalman::Exec hits: " << nHits << std::endl;
            std::cout << "MpdTpcAlignmentKalman::Exec tracks: " << tracks.size() << std::endl;
        }
#endif

        tracksProcess();
    }
    else
        Warning("MpdTpcAlignmentKalman::Exec", "finished with ZERO tracks");

    exec_time += (clock() - tStart);
}


Int_t MpdTpcAlignmentKalman::tracksFind()
{
    FittedTrack tr;
    Int_t sec_hit[nSectors];

    tr.org_mc = Eigen::Vector3d::Zero();
    tr.dir_mc = Eigen::Vector3d::Zero();

    auto track_check_and_store = [&]() -> void
    {
        if (tr.hits.empty())
        {
            tr.clear();
            return;
        }

        Int_t nsec = 0;

        for (unsigned sid = 0; sid < nSectors; ++sid)
            if (sec_hit[sid]) ++nsec;

        if (nsec >= nSectorsPerTrackMin &&  trackParsCompute(tr))
        {
            tracks.push_back(tr);
            tr.reset(); // DO NOT use clear() while storing to vector
        }
        else
            tr.clear();
    };

    if ((!std::strcmp(tracks_type, "laser_mc")) || (!std::strcmp(tracks_type, "laser_emu")))
    {
        Eigen::Vector3d org, dir, hit_pos, dp;

        // collecting hits - laser_mc only
        const double R = 123 - 3. / std::cos(15.*TMath::DegToRad()); // Membrane_outer_holder_R_edge
        const double org_z_offset = 30.; // planes offset, cm
        const double mir_ang = 11.;      // mirror ang, degree
        const double bndl_ang = 8.;      // bundle ang, degree

        double Rc = 1.5;                 // max hit-to-the-track distance, cm
        Rc *= Rc;

        double Rr = 10.;                // min hit-to-the-source distance, cm
        Rr *= Rr;

        auto laser_hits_collect = [&]() -> void
        {
            for (unsigned sid = 0; sid < nSectors; ++sid)
                sec_hit[sid] = 0;

            for (auto it = hits.begin(); it != hits.end(); ++it)
            {
                hit_pos[0] = (*it)->glb.x();
                hit_pos[1] = (*it)->glb.y();
                hit_pos[2] = (*it)->glb.z();

                dp = hit_pos - org;
                auto ds = (dp - (dir.adjoint()*dp)*dir).squaredNorm();
                auto dr = dp.squaredNorm();

                if (ds < Rc && dr > Rr)
                {
                    tr.hits.push_back(*it);
                    ++sec_hit[(*it)->sid];
                }
            }

            tr.org_mc = org;
            tr.dir_mc = dir;

            track_check_and_store();
        };

        for (Int_t pln = -4; pln < 4; ++pln)           // planes
        {
            org[2] = pln < 0 ? pln * org_z_offset : (pln + 1) * org_z_offset;

            for (Int_t beam = 0; beam < 4; ++beam)     // beams
            {
                auto ang = (beam * 90. + bndl_ang) * TMath::DegToRad();
                org[0] = R * std::sin(ang);
                org[1] = R * std::cos(ang);

                for (Int_t refl = -2; refl < 2; ++refl) // reflections
                {
                    auto phi = (360. - 90.*(beam + 1) - bndl_ang + (refl < 0 ? (refl - 1) * mir_ang : (refl + 2) * mir_ang)) * TMath::DegToRad();
                    auto sn  = std::sin(phi);
                    auto cs  = std::cos(phi);

                    dir << cs, sn, 0.;
                    laser_hits_collect();

                    if (0)//(refl > -2)
                    {
                        phi = (360. - 90. * (beam + 1) - bndl_ang + refl * mir_ang) * TMath::DegToRad();
                        sn  = std::sin(phi);
                        cs  = std::cos(phi);

                        dir << cs, sn, 0.;
                        laser_hits_collect();
                    }
                }
            }
        }
    } // laser_mc
    else if (!std::strcmp(tracks_type, "emulation"))
    {
        for (unsigned sid = 0; sid < nSectors; ++sid)
            sec_hit[sid] = 0;

        for (auto it = hits.begin(); it != hits.end(); ++it)
        {
            tr.hits.push_back(*it);
            ++sec_hit[(*it)->sid];
        }

        track_check_and_store();
    } // emulation (1 single track per event)
    else
    {
        // loop over data tracks (external reconstruction, no-Bfield only)
        {
            for (unsigned sid = 0; sid < nSectors; ++sid)
                sec_hit[sid] = 0;

            // assign hits to track, count num of hits per sector segment
            // TODO

            track_check_and_store();
        }
    } // real data

    return tracks.size();
}


bool MpdTpcAlignmentKalman::trackParsCompute(FittedTrack &tr)
{
    bool result = false;

    // Least-squares fit of a line to (x,y,z) data by using distance measurements
    // orthogonal to the proposed line. The line is (P,D) = (origin = mean,
    // direction = the eigenvector corresponding the largest eigenvalue).
    // The error for S = (x0,y0,z0) is (S-P)^T*(I - D*D^T)*(S-P).

    // Compute the mean of the points
    Eigen::Vector3d mean = Eigen::Vector3d::Zero();
    Hit *hit0 = *tr.hits.begin();

    for (auto it = tr.hits.begin(); it != tr.hits.end(); ++it)
    {
        mean += (*it)->glb;

        if ((hit0->time > (*it)->time))// && (hit0->glb.squaredNorm() < (*it)->glb.squaredNorm()))  // find entry point
            hit0 = (*it);
    }

    double invSize = 1./tr.hits.size();
    mean *= invSize;

    if (std::isfinite(mean[0]) && std::isfinite(mean[1]) && std::isfinite(mean[2]))
    {
        // Fit -------------------------------------------------------------

        Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
        Eigen::Vector3d dp;

        for (auto it = tr.hits.begin(); it != tr.hits.end(); ++it)
        {
            dp = (*it)->glb - mean;

            cov(0,0) += dp[0] * dp[0];     cov(0,1) += dp[0] * dp[1];   cov(0,2) += dp[0] * dp[2];
            /*cov(1,0) += dp[1] * dp[0];*/ cov(1,1) += dp[1] * dp[1];   cov(1,2) += dp[1] * dp[2];
            /*cov(2,0) += dp[2] * dp[0];   cov(2,1) += dp[2] * dp[1];*/ cov(2,2) += dp[2] * dp[2];
        }

        cov(0,0) *= invSize; cov(0,1) *= invSize; cov(0,2) *= invSize;
        cov(1,0) = cov(0,1); cov(1,1) *= invSize; cov(1,2) *= invSize;
        cov(2,0) = cov(0,2); cov(2,1) = cov(1,2); cov(2,2) *= invSize;

        // Compute the eigenvector corresponding to the maximum eigenvalue.
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3> > es(cov);
        Eigen::Vector3d::Index maxRow, maxCol;
        es.eigenvalues().maxCoeff(&maxRow, &maxCol);

        tr.dir = es.eigenvectors().col(maxRow);
        tr.org = mean + (tr.dir.adjoint()*(hit0->glb - mean))*tr.dir; // recalculate origin

        //std::cout << " rec org " << tr.org.transpose() << " dir " << tr.dir.transpose() << " hits " << tr.hits.size() << std::endl;

        // Track pars ------------------------------------------------------

        Int_t  side = 0;
        std::vector<Hit*> sector_hits[nSectors];

        for (auto it = tr.hits.begin(); it != tr.hits.end(); ++it)
        {
            auto hit_pos = (*it)->glb;

            tr.side = hit_pos[2] < 0. ? -1 : 1;

            if (side && (side*tr.side < 0))
            {
                Warning("MpdTpcAlignmentKalman::trackFit", "track containes hits from both detector sides; skipping");
                return result;
            }
            else
                side = tr.side;

            auto sid = (*it)->sid;

            if (!tr.cov[sid])
            {
                sector_hits[sid].clear();

                tr.cov[sid] = new Eigen::Matrix3d();
                tr.cov[sid]->setZero();
            }

            sector_hits[sid].push_back(*it);

            auto pos = tr.org + (tr.dir.adjoint()*(hit_pos - tr.org))*tr.dir; // fitted track's point closest to the hit, glob coord
            (*it)->fit = pos;

            auto& tr_cov = (*tr.cov[sid]);
            dp  = hit_pos - pos;
            //std::cout << " dist " << dp.transpose() << " norm " << dp.norm() << " cm " << std::endl;

            tr_cov(0,0) += dp[0] * dp[0];     tr_cov(0,1) += dp[0] * dp[1];   tr_cov(0,2) += dp[0] * dp[2];
            /*tr_cov(1,0) += dp[1] * dp[0];*/ tr_cov(1,1) += dp[1] * dp[1];   tr_cov(1,2) += dp[1] * dp[2];
            /*tr_cov(2,0) += dp[2] * dp[0];   tr_cov(2,1) += dp[2] * dp[1];*/ tr_cov(2,2) += dp[2] * dp[2];
        }

        tr.nsec = 0;

        for (unsigned sid = 0; sid < nSectors; ++sid)
        {
            if (tr.cov[sid])
            {
                auto nshits = sector_hits[sid].size();

                if (nshits < nHitsPerSecMin)
                {
                    delete tr.cov[sid];
                    tr.cov[sid] = nullptr;

                    continue;
                }

                ++tr.nsec;
                invSize = 1./nshits;

                auto& tr_cov = (*tr.cov[sid]);

                tr_cov(0,0) *= invSize;    tr_cov(0,1) *= invSize;    tr_cov(0,2) *= invSize;
                tr_cov(1,0) = tr_cov(0,1); tr_cov(1,1) *= invSize;    tr_cov(1,2) *= invSize;
                tr_cov(2,0) = tr_cov(0,2); tr_cov(2,1) = tr_cov(1,2); tr_cov(2,2) *= invSize;

                /* test
                for (auto tt = sector_hits[sid].begin(); tt != sector_hits[sid].end(); ++tt)
                {
                    auto* hit = (*tt);

                    auto a = hit->glb;
                    auto b = fValue_inv(a, dt);
                    auto c = fValue(b.second, b.first, dt);

                    //auto b = glob2loc(a);
                    //auto c = loc2glob(b.first, b.second);

                    std::cout << (c-a).transpose() << " -- " << (b.first-hit->loc).transpose() << std::endl;
                } // */
#ifdef DEBUG
                hits_per_seg->Fill(nshits);
#endif
            }
        }

        if (tr.nsec >= nSectorsPerTrackMin)
        {
            // sort hits ("silly sort", TODO)
            bool swas_done;

            do
            {
                swas_done =  false;

                for (auto it = tr.hits.begin(); it != tr.hits.end();)
                {
                    auto *hp1 = (*it);

                    if (++it != tr.hits.end())
                    {
                        auto *hp2 = (*it);

                        if (hp1->time > hp2->time)
                        {
                            auto *tp  = hp1;
                            (*(--it)) = hp2;
                            (*(++it)) = tp;

                            swas_done = true;
                        }
                    }
                    else
                        break;
                }
            }
            while (swas_done);

            //for (auto it = tr.hits.begin(); it != tr.hits.end(); ++it)
                //std::cout << (*it)->time << "  " << (*it)->sid << std::endl;

            // check direction sign
            {
                auto *hp1 = *(tr.hits.begin()), *hp2 = *(--tr.hits.end());

                if ((hp2->glb - hp1->glb).dot(tr.dir) < 0)
                {
                    tr.dir *= -1.;
                    //std::cout << " dir inversion " << std::endl;
                }
            }

            tr.chi2 = 0;
            auto nshits = 0;

            for (auto it = tr.hits.begin(); it != tr.hits.end(); ++it)
            {
                auto tr_cov = tr.cov[(*it)->sid];

                if (tr_cov)
                {
                    dp = (*it)->glb - (*it)->fit;
                    tr.chi2 += dp.adjoint()*(*tr_cov)*dp;
                    ++nshits;
                }
                //else
                    //it = tr.hits.erase(it);
            }

            tr.chi2 /= nshits;

            result = true;
        }

#ifdef DEBUG
        hits_per_trk->Fill(tr.hits.size());
        segs_per_trk->Fill(tr.nsec);
#endif
    }

    return result;
}


void MpdTpcAlignmentKalman::tracksProcess()
{
    Eigen::Matrix<double, nsums, 1>     sn; // predicted
    Eigen::Matrix<double, nsums, nsums> Cn; // predicted

    Eigen::Matrix<double, ndims, 1> vel;
    Eigen::Matrix<double, ndims, 1> res;
    Eigen::Matrix<double, nsums, ndims> K;
    Eigen::Matrix<double, ndims, ndims> R, Rinv;

    bool invertible, sector_active;
    unsigned sid_n, sid_p;
    Hit *hit_n, *hit_p;

    for (auto tr = tracks.begin(); tr != tracks.end(); ++tr)
    {
        auto it = tr->hits.begin();
        vel = (tr->dir.block(0,0,ndims,1).normalized().array()*c_vel).matrix();

        hit_p = *it;
        sid_p = hit_p->sid;

        sector_active = sector_state_init(sid_p, vel, hit_p->fit.block(0,0,ndims,1), &(*tr));

        for (++it; it != tr->hits.end(); ++it)
        {
            hit_n = *it;
            sid_n = hit_n->sid;

            // sector swap
            if (sid_n != sid_p)
            {
                // save current alignment state for the previous sector
                if (sector_active)
                    sector_state_revise(sid_p);

                sector_active = sector_state_init(sid_n, vel, hit_n->fit, &(*tr));

                hit_p = hit_n;
                sid_p = sid_n;

                continue;
            }

            // process state vector
            if (sector_active)
            {
                // prediction
                auto dt = hit_n->time - hit_p->time;
                A.block(0,ndims,ndims,ndims) = (Ind.array()*dt).matrix();

                //sn = A*sp;
                sn = sp; sn.block(0,0,ndims,1) = hit_n->fit.block(0,0,ndims,1);
                Cn = A*Cp*A.transpose() + Q;

                R = V + H[sid_n]*Cn*H[sid_n].transpose();
                R.computeInverseWithCheck(Rinv, invertible);

                // correction
                if (invertible)
                {
                    res = hit_n->raw.block(0,0,ndims,1) - H[sid_n]*sn;
                    K   = Cp*H[sid_n].transpose()*Rinv;
                    /*
                    std::cout << "+ dlprd,cm: " << (sn.block(0,0,ndims,1) - sp.block(0,0,ndims,1)).norm() << " vdt: "  << c_vel*dt << std::endl;
                    std::cout << "  line dir: " << (hit_n->nrm.block(0,0,ndims,1) - hit_p->nrm.block(0,0,ndims,1)).normalized().transpose() <<
                                 " v dir: "  << (sp.block(ndims,0,ndims,1)*dt).normalized().transpose() << std::endl;
                    std::cout << std::endl << K << " K" << std::endl << std::endl;
                    std::cout << "  mhit: " << hit_n->glb.block(0,0,ndims,1).transpose() << " model: "  << (H[sid_n]*sn).transpose() << std::endl;
                    std::cout << "  res:  " << res.transpose() << " norm: "  << res.norm() << std::endl << K*res << " K*res" << std::endl;
                    std::cout << "  sold: " << sp.transpose() << std::endl << "  sprd: " << sn.transpose() << std::endl; */

                    sp = sn + K*res;
                    Cp = (Ins - K*H[sid_n])*Cn;

                    //std::cout << "  snew: " << sp.transpose() << std::endl;

                    R = V + H[sid_n]*Cp*H[sid_n].transpose();
                    R.computeInverseWithCheck(Rinv, invertible);

                    if (invertible)
                        chi2[sid_n] += res.transpose()*Rinv*res;
                }
            }

            hit_p = hit_n;
            sid_p = sid_n;
        }

        if (sector_active)
            sector_state_revise(sid_p);

        ++tracks_processed;
    }

    ++events_processed;

#ifdef DEBUG
    event_log_flush();
#endif
}


void MpdTpcAlignmentKalman::matrixHInit(const unsigned int& sid)
{
    // measurement:
    // 2loc -> misalign -> 2glob -> hit coord

    auto k = dPhi*sid;
    auto cs = std::cos(k);
    auto sn = std::sin(k);

    Eigen::Matrix<double, ndims, ndims> h1, h2, h3;
    Eigen::Matrix<double, nsums, nsums> T1, T2;
    Eigen::Matrix<double, ndims, nsums> T3;

    h1 << cs, sn, 0, -sn, cs, 0, 0, 0, 1;
    h2 = (Ind.array()*(-1.)).matrix();
    h3 << cs, -sn, 0, sn, cs, 0, 0, 0, 1;

    T1.setZero().block(0,0,ndims,ndims) = h1;
    T1.block(ndims,ndims,ndims,ndims)   = Ind;
    T1.block(nipar,nipar,npars,npars)   = Inp;

    T2.setZero().block(0,0,ndims,ndims) = Ind;
    T2.block(0,nipar,ndims,ndims)       = h2;
    T2.block(ndims,ndims,ndims,ndims)   = Ind;
    T2.block(nipar,nipar,npars,npars)   = Inp;

    T3.setZero().block(0,0,ndims,ndims) = h3;

    H[sid] = T3*T2*T1;
}


bool MpdTpcAlignmentKalman::sector_state_init(const unsigned int &sid, const Eigen::Matrix<double, ndims, 1> &vel, const Eigen::Matrix<double, ndims, 1> &pos, const FittedTrack *tr)
{
    // define if current sector is to be processed
    bool sector_active = sector_status_(sid) && tr->cov[sid];

    if (sector_active)
    {
        weak_modes_fix(s);

        sp.block(0,0,ndims,1)     = pos;
        sp.block(ndims,0,ndims,1) = vel;
        sp.block(nipar,0,npars,1) = s[sid];

        Cp.setZero().diagonal().setConstant(Cini);
        Cp.block(nipar,nipar,npars,npars) = C[sid];

        //std::cout << "  sector " << sid << std::endl  << std::endl << H[sid] << " H" << std::endl << std::endl;

        ann[sid] = (track_segments[sid] < nAnnSecTracks ? std::pow(ann0, static_cast<double>(nAnnSecTracks-track_segments[sid])/static_cast<double>(nAnnSecTracks)) : 1.);
        auto annCoeff = Ind.array()*ann[sid];

        V = ((*tr->cov[sid]).array()*annCoeff).matrix();
        //V = (*tr->cov[sid]);
        //V.setZero().diagonal().setConstant(1.E0);

        //std::cout << V << " V" << std::endl << std::endl;
    }

    return sector_active;
}


void MpdTpcAlignmentKalman::sector_state_revise(const unsigned int &sid)
{
    ++track_segments[sid];

    s[sid] = sp.block(nipar,0,npars,1);
    C[sid] = Cp.block(nipar,nipar,npars,npars);

    pstats[sid].add(s[sid]);
}


bool MpdTpcAlignmentKalman::sector_status_(const unsigned int &sid)
{
    /* check if alignment of the sector is over
    if (!is_finished[sid])
    {
        bool is_over = true;

        for (auto i = 0; i < npars; ++i)
            is_over = (is_over && (C[sid](i,i) < Cmin));

        if (is_over) is_finished[sid] = true;
    } */

    return is_active[sid] && !is_finished[sid];
}


void MpdTpcAlignmentKalman::weak_modes_fix(Eigen::Matrix<double, npars, 1> *pars)
{
    //return;

    Eigen::Vector3d cmean, s3;
    double pmean;
    int nsids;

    for (unsigned cam = 0; cam < 2; ++cam)
    {
        cmean.setZero();
        pmean = 0;
        nsids = 0;

        for (unsigned sid = cam*nSectors_half; sid < (cam+1)*nSectors_half; ++sid)
        {
            if (is_active[sid] && !is_finished[sid])
            {
                s3.setZero().block(0,0,ndims,1) = pars[sid].block(0,0,ndims,1);

                cmean += loc2glob(s3, sid);
                //pmean += pars[sid](p_dph,1);

                ++nsids;
            }
        }

        if (!nsids) continue;
        double norm = 1./nsids;

        cmean *= norm;
        pmean *= norm;

        for (unsigned sid = cam*nSectors_half; sid < (cam+1)*nSectors_half; ++sid)
        {
            if (is_active[sid] && !is_finished[sid])
            {
                s3.setZero().block(0,0,ndims,1) = pars[sid].block(0,0,ndims,1);

                pars[sid].block(0,0,ndims,1) = glob2loc((loc2glob(s3, sid) - cmean), sid).first.block(0,0,ndims,1);
                //pars[sid](p_dph,1) -= pmean;
            }
        }
    }

}


MpdTpcAlignmentKalman::Hit* MpdTpcAlignmentKalman::hit_(const MpdTpcHit* tpcHit) const
{
    Hit* hit = new Hit();
    TVector3 globPos;

    hit->time = tpcHit->GetTimeStamp();
    tpcHit->Position(globPos);

    hit->raw << globPos.X(), globPos.Y(), globPos.Z();
    auto locPos = glob2loc(hit->raw);

    hit->loc = locPos.first;
    hit->sid = locPos.second;
    hit->glb = /*track_segments[hit->sid] < nAnnSecTracks ? hit->raw :*/ align(locPos.second, locPos.first); // apply currently known alignment for hits to be used in fit

    //std::cout << "sector: " << (hit)->sid << " glob: " << (hit)->glb.transpose() << " loc: " << (hit)->loc.transpose() << std::endl;
    //std::cout << "sector: " << sid << " glob X: " << globPos.X() << " Y: " << globPos.Y() << " Z: " << globPos.Z() << " loc X: " << locPos.first.X() << " Y: " << locPos.first.Y() << " Z: " << locPos.first.Z() << std::endl;

    return hit;
}


unsigned MpdTpcAlignmentKalman::sid_(const Eigen::Vector3d& globPos) const
{
    auto phi = std::atan2(globPos(1,0), globPos(0,0)) - M_PI_2;
    while (phi < 0) phi += 2.*M_PI;

    auto sid = (int)(phi/dPhi + 0.5);
    if (sid == nSectors_half) sid = 0;

    return static_cast<unsigned>(globPos(2,0) > 0.0 ? sid : sid + nSectors_half);
}


std::pair<Eigen::Vector3d, unsigned> MpdTpcAlignmentKalman::glob2loc(const Eigen::Vector3d& globPos, const unsigned sid2use) const
{
    unsigned sid = sid2use;
    Eigen::Vector3d locPos;

    if (sid2use >= nSectors)
        sid = sid_(globPos);

    auto k = dPhi*sid;
    auto cs = std::cos(k);
    auto sn = std::sin(k);

    auto x = globPos(0,0);
    auto y = globPos(1,0);

    locPos(0,0) = x*cs + y*sn;
    locPos(1,0) = y*cs - x*sn;// - lowEdgeY;
    locPos(2,0) = globPos(2,0);// + (sid < nSectors_half ? -driftLength : driftLength);

    return std::pair<Eigen::Vector3d, unsigned>(locPos, sid);
}


Eigen::Vector3d MpdTpcAlignmentKalman::loc2glob(const Eigen::Vector3d& locPos, const unsigned& sid) const
{
    Eigen::Vector3d globPos;
    globPos(2,0) = locPos(2,0);// + (sid < nSectors_half ? driftLength : -driftLength);

    auto k = dPhi*sid;
    auto cs = std::cos(-k);
    auto sn = std::sin(-k);

    auto x = locPos(0,0);
    auto y = locPos(1,0);//+ lowEdgeY;

    globPos(0,0) = x*cs + y*sn;
    globPos(1,0) = y*cs - x*sn;

    return globPos;
}


Eigen::Vector3d MpdTpcAlignmentKalman::align(const unsigned& sid, const Eigen::Vector3d& arg, const Eigen::Matrix<double, npars, 1>* a) const
{
    Eigen::Vector3d val;
    const Eigen::Matrix<double, npars, 1>* pars = (a ? a : s);

    //auto p = pars[sid](p_dp,0);
    //auto sn = std::sin(p);
    //auto cs = std::cos(p);

    auto x = arg(0,0);
    auto y = arg(1,0);
    auto z = arg(2,0);

    //val(0,0) = x*cs+y*sn;
    //val(1,0) = y*cs-x*sn;

    val(0,0) = x + pars[sid](p_dx,0);
    val(1,0) = y + pars[sid](p_dy,0);
    val(2,0) = z + pars[sid](p_dz,0);

    return loc2glob(val, sid);
}


std::pair<Eigen::Vector3d, unsigned> MpdTpcAlignmentKalman::misalign(const Eigen::Vector3d &arg, const Eigen::Matrix<double, npars, 1>* a) const
{
    auto  loc = glob2loc(arg);
    auto  sid = loc.second;
    auto& val = loc.first;

    const Eigen::Matrix<double, npars, 1>* pars = (a ? a : s);

    //auto p = pars[sid](p_dp,0);
    //auto sn = std::sin(-p);
    //auto cs = std::cos(-p);

    auto x = val(0,0) - pars[sid](p_dx,0);
    auto y = val(1,0) - pars[sid](p_dy,0);
    auto z = val(2,0) - pars[sid](p_dz,0);

    //val(0,0) = x*cs+y*sn;
    //val(1,0) = y*cs-x*sn;
    //val(2,0) = z;

    val(0,0) = x;  // no dph
    val(1,0) = y;  //
    val(2,0) = z;  //

    return loc;
}


void MpdTpcAlignmentKalman::FinishEvent()
{
    if (hits.empty())
        return;

    // cleanup
    for (auto it = hits.begin(); it != hits.end(); ++it)
        if (*it) { delete *it; *it = nullptr; }

    for (auto tr = tracks.begin(); tr != tracks.end(); ++tr)
        tr->clear();
}


void MpdTpcAlignmentKalman::Finish()
{
    auto time_tot = ((float) exec_time) / CLOCKS_PER_SEC;
    std::cout << "MpdTpcAlignmentKalman execution time: " << time_tot << " s, per event: " << time_tot / events_processed << " s" << std::endl;

#ifdef DEBUG
    out_file->cd();
    out_tree->Write();

    hits_per_trk->Write();
    hits_per_seg->Write();
    segs_per_trk->Write();

    for (unsigned sid = 0; sid < nSectors; ++sid)
    {
        hits_per_sec[sid]->Write();
        segs_per_sec[sid]->Write();
    }

    out_file->Write();
    //out_file->WriteTObject(out_tree);
    out_file->Close();

    //if (reader)
    //{
    //    delete reader; reader = nullptr;
    //}

    if (rand_gen)
    {
        delete rand_gen; rand_gen = nullptr;
    }
#endif
}

#ifdef DEBUG
double MpdTpcAlignmentKalman::rnd_flat(const double min, const double max) const
{
    std::uniform_real_distribution<double> range(min, max);
    return range(*rand_gen);
}


double MpdTpcAlignmentKalman::rnd_norm(const double mean, const double sigma) const
{
    if (sigma > 0)
    {
        std::normal_distribution<double> range(mean, sigma);
        return range(*rand_gen);
    }

    return mean;
}


double MpdTpcAlignmentKalman::rnd_exp(const double mean) const
{
    std::exponential_distribution<double> range(mean);
    return range(*rand_gen);
}


void MpdTpcAlignmentKalman::misalignment_commit(unsigned int nHits)
{
    for (unsigned i = 0; i < nHits; ++i)
    {
        MpdTpcHit * hit = (MpdTpcHit*)fHits->UncheckedAt(i);

        TVector3 globPos;
        hit->Position(globPos);

        Eigen::Vector3d glob;
        glob << globPos.X(), globPos.Y(), globPos.Z();

        auto loc = misalign(glob, st);

        glob = loc2glob(loc.first, loc.second);
        globPos.SetXYZ(glob(0,0), glob(1,0), glob(2,0));

        TVector3 locPos;
        fSecGeo->Global2Local(globPos, locPos);

        hit->SetPosition(globPos);
        hit->SetLocalPosition(locPos);
    }
}


void MpdTpcAlignmentKalman::event_log_flush()
{
    //return;

    Int_t nshit[nSectors], nsseg[nSectors];

    for (unsigned sid = 0; sid < nSectors; ++sid)
        nshit[sid] = nsseg[sid] = 0;

    for (auto tr = tracks.begin(); tr != tracks.end(); ++tr)
    {
        for (unsigned sid = 0; sid < nSectors; ++sid)
        {
            if (tr->cov[sid])
            {
                ++nsseg[sid];

                for (auto it = tr->hits.begin(); it != tr->hits.end(); ++it)
                {
                    if (sid == (*it)->sid)
                        ++nshit[sid];
                }
            }
        }
    }

    Eigen::Matrix<double, npars, 1> mean, disp;

    if ( !(events_processed%printout_period))
    {
        std::cout << "Event #" << events_processed << ": sec0 ann " << ann[0] <<  " Cxx " <<  C[0](p_dx,p_dx) <<  " Cyy " <<  C[0](p_dy,p_dy) <<  " Czz " <<  C[0](p_dz,p_dz);
        std::cout << std::endl << "                   diff-diff_prev                diff/true                     true                  current                       sigma    tracks" << std::endl;
    }

    for (unsigned sid = 0; sid < nSectors; ++sid)
    {
        bool stats_ready = false; //track_segments[sid] > 1 ? true : false;
        pstats[sid].stats_compute(mean, disp);

        ost.t  [sid] = events_processed;
        ost.ann[sid] = ann[sid];
        ost.nsg[sid] = nsseg[sid];
        ost.nht[sid] = nshit[sid];
        ost.dx [sid] = (stats_ready ? mean(p_dx,0) : s[sid](p_dx,0)) - st[sid](p_dx,0);
        ost.dy [sid] = (stats_ready ? mean(p_dy,0) : s[sid](p_dy,0)) - st[sid](p_dy,0);
        ost.dz [sid] = (stats_ready ? mean(p_dz,0) : s[sid](p_dz,0)) - st[sid](p_dz,0);
        ost.dp [sid] = 0.;//(stats_ready ? mean(p_dp,0) : s[sid](p_dp,0)) - st[sid](p_dp,0);
        ost.de [sid] = 0.;//(stats_ready ? mean(p_de,0) : s[sid](p_de,0)) - st[sid](p_de,0);
        ost.dk [sid] = 0.;//(stats_ready ? mean(p_dk,0) : s[sid](p_dk,0)) - st[sid](p_dk,0);
        ost.Cx [sid] = C[sid](p_dx,p_dx);
        ost.Cy [sid] = C[sid](p_dy,p_dy);
        ost.Cz [sid] = C[sid](p_dz,p_dz);
        ost.Cp [sid] = 0.;//C[sid](p_dp,p_dp);
        ost.Ce [sid] = 0.;//C[sid](p_de,p_de);
        ost.Ck [sid] = 0.;//C[sid](p_dk,p_dk);
        //ost.Vx [sid] = V(s_x,s_x);
        //ost.Vy [sid] = V(s_y,s_y);
        //ost.Vz [sid] = V(s_z,s_z);

        auto diffx = std::abs(ost.dx[sid]);
        auto diffy = std::abs(ost.dy[sid]);
        auto diffz = std::abs(ost.dz[sid]);
        auto qx = diffx/std::abs(st[sid](p_dx,0));
        auto qy = diffy/std::abs(st[sid](p_dy,0));
        auto qz = diffz/std::abs(st[sid](p_dz,0));

        if (is_active[sid] &&  !(events_processed%printout_period))
            std::printf("sect %02d:     % .0e % .0e % .0e     % .3f % .3f % .3f     % .3f % .3f % .3f     % .3f % .3f % .3f     %.1e %.1e %.1e     %5d\n", sid,
                        (diffx-s0[sid](p_dx,0)), (diffy-s0[sid](p_dy,0)), (diffz-s0[sid](p_dz,0)),
                                         qx,                      qy,                      qz,
                               st[sid](p_dx,0),         st[sid](p_dy,0),         st[sid](p_dz,0),
                                s[sid](p_dx,0),          s[sid](p_dy,0),          s[sid](p_dz,0),
                        std::sqrt(disp(p_dx,0)), std::sqrt(disp(p_dy,0)), std::sqrt(disp(p_dz,0)), track_segments[sid]  );

        s0[sid](p_dx,0) = diffx;
        s0[sid](p_dy,0) = diffy;
        s0[sid](p_dz,0) = diffz;

        if (events_processed)
        {
            hits_per_sec[sid]->Fill(nshit[sid]);
            segs_per_sec[sid]->Fill(nsseg[sid]);
        }
    }

    out_tree->Fill();
}


void MpdTpcAlignmentKalman::matrix_print(unsigned sid)
{
    std::cout << "--------- Sect " << sid << std::endl; // << " ds " << (s[sid] - st[sid]).transpose() << std::endl;
    std::cout << " curr " << (s[sid]).transpose() << std::endl;
    //std::cout << " sprd " << sprd.transpose() << std::endl << std::endl;
    //std::cout << Cprd  << " Cprd" << std::endl << std::endl;
    //std::cout << K  << " K" << std::endl << std::endl;
    std::cout << A  << " A" << std::endl << std::endl;
    std::cout << H[sid]  << " H" << std::endl << std::endl;
    std::cout << C[sid]  << " C" << std::endl << std::endl;
    //std::cout << Q  << " D" << std::endl << std::endl;
    //std::cout << V  << " V" << std::endl << std::endl;
    std::cout << "----------------------" << std::endl;
}


/*
//std::pair<TVector3, unsigned> MpdTpcAlignmentKalman::glob2loc(const TVector3& globPos) const
std::pair<Eigen::Vector3d, unsigned> MpdTpcAlignmentKalman::glob2loc(const TVector3& globPos) const
{
    Eigen::Vector3d glob;
    glob << globPos.X(), globPos.Y(), globPos.Z();

    auto loc = glob2loc(glob);
    //return std::pair<TVector3, unsigned>(TVector3(loc.first(0,0), loc.first(1,0), loc.first(2,0)), loc.second);
    return loc;
}


//TVector3 MpdTpcAlignmentKalman::loc2glob(const TVector3& locPos, const unsigned& sid) const
Eigen::Vector3d MpdTpcAlignmentKalman::loc2glob(const TVector3& locPos, const unsigned& sid) const
{
    Eigen::Vector3d loc;
    loc << locPos.X(), locPos.Y(), locPos.Z();

    auto glob = loc2glob(loc, sid);
    //return TVector3(glob(0,0), glob(1,0), glob(2,0));
    return glob;
}


void MpdTpcAlignmentKalman::hit2vector(const Hit* hit, TVector3& globXYZ, TVector3& locXYZ) const
{
    globXYZ.SetX(hit->glb(0));
    globXYZ.SetY(hit->glb(1));
    globXYZ.SetZ(hit->glb(2));

    fSecGeo->Global2Local(globXYZ, locXYZ); // local CS defined in MpdTpcSectorGeo
} */
#endif

ClassImp(MpdTpcAlignmentKalman)
