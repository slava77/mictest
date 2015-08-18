#ifndef CandCloner_h
#define CandCloner_h

#include "SideThread.h"

#include "MkFitter.h"

#include <vector>

//#define CC_TIME_LAYER
//#define CC_TIME_ETA

typedef std::pair<int, int> CandClonerWork_t;


class CandCloner : public SideThread<CandClonerWork_t>
{
public:
  // Maximum number of seeds processed in one call to ProcessSeedRange()
  static const int s_max_seed_range = 4;

private:
  // Temporaries in ProcessSeedRange(), re-sized/served  in constructor.

  // Size of this guy is s_max_seed_range * max_cand
  std::vector<std::pair<int,MkFitter::IdxChi2List> > t_seed_newcand_idx;

  // Size of this one is s_max_seed_range
  std::vector<std::vector<Track> > t_cands_for_next_lay;

public:
  CandCloner(int cpuid=-1, int cpuid_st=-1)
  {
    m_fitter = new (_mm_malloc(sizeof(MkFitter), 64)) MkFitter(0);

    t_seed_newcand_idx.reserve(s_max_seed_range * Config::maxCand);

    t_cands_for_next_lay.resize(s_max_seed_range);
    for (int iseed = 0; iseed < s_max_seed_range; ++iseed)
    {
      t_cands_for_next_lay[iseed].reserve(Config::maxCand);
    }

#ifndef CLONE_ENGINE_SINGLE_THREAD
    SpawnSideThread(cpuid, cpuid_st);
#endif
  }

  ~CandCloner()
  {
    // printf("CandCloner::~CandCloner will try to join the side thread now ...\n");

#ifndef CLONE_ENGINE_SINGLE_THREAD
    JoinSideThread();
#endif

    _mm_free(m_fitter);
  }

  void begin_eta_bin(EtaBinOfCombCandidates * eb_o_ccs, int start_seed, int n_seeds)
  {
    mp_etabin_of_comb_candidates = eb_o_ccs;
    m_start_seed = start_seed;
    m_n_seeds    = n_seeds;
    m_hits_to_add.resize(n_seeds);

    // XXX Should resize vectors in m_hits_to_add to whatever makes sense ???
    // for (int i = 0; i < n_seeds; ++i)
    // {
    //   m_hits_to_add[i].reserve(20 * Config::maxCand);
    // }

#ifdef CC_TIME_ETA
    printf("CandCloner::begin_eta_bin\n");
    t_eta = dtime();
#endif
  }

  void begin_layer(BunchOfHits *b_o_hs, int lay)
  {
    mp_bunch_of_hits = b_o_hs;

    m_layer = lay;

    m_fitter->SetNhits(m_layer);//here again assuming one hit per layer

    m_idx_max      = 0;
    m_idx_max_prev = 0;

#ifdef CC_TIME_LAYER
    t_lay = dtime();
#endif
  }

  void begin_iteration()
  {
    // Do nothing, "secondary" state vars updated when work completed/assigned.
  }

  void add_cand(int idx, const MkFitter::IdxChi2List& cand_info)
  {
    m_hits_to_add[idx].push_back(cand_info);

    m_idx_max = std::max(m_idx_max, idx);
  }

  void end_iteration()
  {
    int proc_n = m_idx_max - m_idx_max_prev;

    // printf("CandCloner::end_iteration process %d, max_prev=%d, max=%d\n", proc_n, m_idx_max_prev, m_idx_max);

    if (proc_n >= s_max_seed_range)
    {
      // Round to multiple of s_max_seed_range.
      signal_work_to_st((m_idx_max / s_max_seed_range) * s_max_seed_range);
    }
  }

  void end_layer()
  {
    if (m_idx_max + 1 > m_idx_max_prev)
    {
      signal_work_to_st(m_idx_max + 1);
    }

#ifndef CLONE_ENGINE_SINGLE_THREAD
    WaitForSideThreadToFinish();
#endif

    for (int i = 0; i < m_n_seeds; ++i)
    {
      m_hits_to_add[i].clear();
    }

#ifdef CC_TIME_LAYER
    t_lay = dtime() - t_lay;
    printf("CandCloner::end_layer %d -- t_lay=%8.6f\n", m_layer, t_lay);
#endif
  }

  void end_eta_bin()
  {
#ifdef CC_TIME_ETA
    t_eta = dtime() - t_eta;
    printf("CandCloner::end_eta_bin t_eta=%8.6f\n", t_eta);
#endif
  }

  void signal_work_to_st(int idx)
  {
    // printf("CandCloner::signal_work_to_st assigning work up to seed %d\n", idx);

#ifndef CLONE_ENGINE_SINGLE_THREAD
    QueueWork(std::make_pair(m_idx_max_prev, idx));
#else
    DoWorkInSideThread(std::make_pair(m_idx_max_prev, idx));
#endif

    m_idx_max_prev = idx;
  }

  // ----------------------------------------------------------------

  void CutAndPastaFromBuildTestMPlex();

  void ProcessSeedRange(int is_beg, int is_end);

  // virtual
  void DoWorkInSideThread(CandClonerWork_t work);

  // ----------------------------------------------------------------

  // eventually, protected or private

  int  m_idx_max, m_idx_max_prev;
  std::vector<std::vector<MkFitter::IdxChi2List>> m_hits_to_add;

  BunchOfHits            *mp_bunch_of_hits;
  EtaBinOfCombCandidates *mp_etabin_of_comb_candidates;

  MkFitter               *m_fitter;

#if defined(CC_TIME_ETA) or defined(CC_TIME_LAYER)
  double    t_eta, t_lay;
#endif

  int       m_start_seed, m_n_seeds;
  int       m_layer;
};

#endif
