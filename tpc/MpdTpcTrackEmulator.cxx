#include "MpdTpcTrackEmulator.h"

MpdTpcTrackEmulator::MpdTpcTrackEmulator(Option_t *emu_tracks_type)
    : FairTask("TPC Track Emulator"), rand_gen(nullptr), fSecGeo(nullptr)
{
    tracks_type = new char[255];
    std::strcpy(tracks_type, emu_tracks_type);

    stddev = mfpath = hitres = 0;
    ntrk_per_ev = 1;
    velocity = 29.9792458; // cm/ns
}


InitStatus MpdTpcTrackEmulator::Init()
{
    if (!(mfpath*stddev*hitres))
    {
        Error("MpdTpcTrackEmulator::Init", "mean free path / scattering standart deviation / hits resolution HAVE TO BE SET explicitly");
        return kFATAL;
    }

    FairRootManager* frman = FairRootManager::Instance();

    if (!frman)
    {
        Error("MpdTpcTrackEmulator::Init", "RootManager not found!");
        return kFATAL;
    }

    fHits = (TClonesArray*) frman->GetObject("TpcRecPoint");

    if (fHits)
    {
        Error("MpdTpcTrackEmulator::Init", "TPC Hits array present");
        return kFATAL;
    }
    else
    {
        fHits = new TClonesArray("MpdTpcHit");
        frman->Register("TpcRecPoint", "TpcAlign", fHits, kTRUE);
    }

    fSecGeo = MpdTpcSectorGeo::Instance();

    if (!fSecGeo)
    {
        Error("MpdTpcTrackEmulator::Init", "TPC Sector Geo init error");
        return kFATAL;
    }

    std::random_device device;
    rand_gen = new std::mt19937(device());

    return kSUCCESS;
}


void MpdTpcTrackEmulator::Exec(Option_t* opt)
{
    hits.clear();

    vertexes_emulate();

    auto nHits = hits.size();
    //std::cout << "MpdTpcTrackEmulator::Exec hits: " << nHits << std::endl;

    if (nHits)
    {
        fHits->Delete();

        for (unsigned i = 0; i < nHits; ++i)
        {
            // reject hits close to the sectors borders
            auto phi = std::atan2(hits[i].first.Y(), hits[i].first.X()) - M_PI_2;
            while (phi < 0) phi += 2.*M_PI;

            auto sid = phi/dPhi + 0.5;
            auto val = sid - (int)(sid);

            if (val < 0.08 || val > 0.92)
                continue;

            // save
            TVector3 locPos;
            fSecGeo->Global2Local(hits[i].first, locPos);

            Int_t iHit = fHits->GetEntriesFast();
            MpdTpcHit* hit = new ((*fHits)[iHit]) MpdTpcHit();

            hit->SetPosition(hits[i].first);
            hit->SetTimeStamp(hits[i].second);
            hit->SetLocalPosition(locPos);
        }
    }
}


void MpdTpcTrackEmulator::vertexes_emulate()
{
    // tracks origin and direction
    // generate multiple vertexes here if needed

    Eigen::Vector3d org, dir;

    if (!std::strcmp(tracks_type, "laser"))
    {
        const double R = R_max - 3. / std::cos(15.*TMath::DegToRad()); // Membrane_outer_holder_R_edge
        const double org_z_offset = 30.; // planes offset, cm
        const double mir_ang = 11.;      // mirror ang, degree
        const double bndl_ang = 8.;      // bundle ang, degree

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
                    track_pass(org, dir);

                    if (0)//(refl > -2)
                    {
                        phi = (360. - 90. * (beam + 1) - bndl_ang + refl * mir_ang) * TMath::DegToRad();
                        sn  = std::sin(phi);
                        cs  = std::cos(phi);

                        dir << cs, sn, 0.;
                        track_pass(org, dir);
                    }
                }
            }
        }
    } // laser
    else if (!std::strcmp(tracks_type, "muons"))
    {
        for (Int_t trk = 0; trk < ntrk_per_ev; ++trk)
        {
            unsigned nHits = 0;

            do
            {
                // 24x4 beams
                auto f = 2, num = 12*f;
                auto k = (int) (rnd_flat(0,num));
                auto ang = 2.*M_PI/num;

                auto phi = ang*k;
                auto ksi = rnd_flat(0,1) > 0.5 ?
                            ang*(rnd_flat(0,1) > 0.5 ? k+f*0.25 : k-f*0.25) :
                            ang*(rnd_flat(0,1) > 0.5 ? k+f*1.50 : k-f*1.50);

                org << R_max*std::cos(phi), R_max*std::sin(phi), (rnd_flat(0,1) > 0.5 ? 80. : -80.);
                dir << -std::cos(ksi), -std::sin(ksi), 0.;

                nHits = track_pass(org, dir);
                //std::cout << " emu org " << org.transpose() << " dir " << dir.transpose() << " hits " << nHits << std::endl;
            }
            while (!nHits);
        }
    } // muons
    else
    {
        for (Int_t trk = 0; trk < ntrk_per_ev; ++trk)
        {
            unsigned nHits = 0;

            do
            {
                // 24x4 beams
                auto f = 2, num = 12*f;
                auto k = (int) (rnd_flat(0,num));
                auto ang = 2.*M_PI/num;

                auto phi = ang*k;
                auto ksi = rnd_flat(0,1) > 0.5 ?
                            ang*(rnd_flat(0,1) > 0.5 ? k+f*0.25 : k-f*0.25) :
                            ang*(rnd_flat(0,1) > 0.5 ? k+f*1.50 : k-f*1.50);

                org << R_max*std::cos(phi), R_max*std::sin(phi), (rnd_flat(0,1) > 0.5 ? 80. : -80.);
                dir << -std::cos(ksi), -std::sin(ksi), 0.;

                /* cross
                if (rnd_flat(0,1) > 0.5)
                {
                    auto R = rnd_flat(R_min, R_max);
                    auto ksi = rnd_flat(0., 2*M_PI);

                    org << 0, R, (rnd_flat(0,1) > 0.5 ? 80. : -80.);
                    dir << std::cos(ksi), std::sin(ksi), 0.;
                }
                else
                {
                    auto R_med = 0.5*(R_max+R_min);
                    auto phi = rnd_flat(M_PI/2 - M_PI/12., M_PI/2 + M_PI/12.);
                    auto ksi = rnd_flat(0., 2*M_PI);

                    org << R_med*std::cos(phi), R_med*std::sin(phi), (rnd_flat(0,1) > 0.5 ? 80. : -80.);
                    dir << std::cos(ksi), std::sin(ksi), 0.;
                } */

                nHits = track_pass(org, dir);
                //std::cout << " emu org " << org.transpose() << " dir " << dir.transpose() << " hits " << nHits << std::endl;
            }
            while (!nHits);
        }
    }  // models
}


unsigned MpdTpcTrackEmulator::track_pass(const Eigen::Vector3d &org, const Eigen::Vector3d &dir)
{
    // track hits emulation

    double hR2 = -1.;
    unsigned nHits = 0;
    Eigen::Vector3d pos = org, hit, dl, ds;

    do
    {
        if (hR2 > R2_min && hR2 < R2_max && fabs(hit(2,0)) < L_max) // register hits within tpc sensitive part only
        {
            // save
            Double_t time = (pos-org).norm()/velocity; // nanoseconds
            hits.push_back(std::pair<TVector3, Double_t>(TVector3(hit(0,0), hit(1,0), hit(2,0)), time));
            ++nHits;
        }

        dl = rnd_exp(mfpath)*dir;

        while (dl.norm() < hitres)
            dl  = rnd_exp(mfpath)*dir;

        pos += dl; // track end position update

        hit = Eigen::Vector3d(rnd_flat(-1, 1), rnd_flat(-1, 1), rnd_flat(-1, 1)); // random vector
        ds  = rnd_norm(0, stddev)*(hit - hit.dot(dir)*dir).normalized();   // perpendicular to dir

        hit = pos + dl + ds;                    // hit position
        hR2 = hit.block(0,0,2,1).squaredNorm(); // hit XY plane r^2
    }
    // widespread outside sources, transparent tpc core
    while (((hR2 < 50.*R2_max && !nHits) || hR2 < R2_max)
           && ((fabs(hit(2,0)) < 50.*L_max && !nHits) || fabs(hit(2,0)) < L_max));
           //&& hR2 > R2_min);

    return nHits;
}


double MpdTpcTrackEmulator::rnd_flat(const double min, const double max) const
{
    std::uniform_real_distribution<double> range(min, max);
    return range(*rand_gen);
}


double MpdTpcTrackEmulator::rnd_norm(const double mean, const double sigma) const
{
    if (sigma > 0)
    {
        std::normal_distribution<double> range(mean, sigma);
        return range(*rand_gen);
    }

    return mean;
}


double MpdTpcTrackEmulator::rnd_exp(const double mean) const
{
    std::exponential_distribution<double> range(mean);
    return range(*rand_gen);
}


ClassImp(MpdTpcTrackEmulator)
