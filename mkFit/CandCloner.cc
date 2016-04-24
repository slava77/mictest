#include "CandCloner.h"

namespace
{
bool sortCandListByHitsChi2(const MkFitter::IdxChi2List& cand1,
                            const MkFitter::IdxChi2List& cand2)
{
  if (cand1.nhits == cand2.nhits) return cand1.chi2 < cand2.chi2;
  return cand1.nhits > cand2.nhits;
}
}

//==============================================================================

void CandCloner::ProcessSeedRange(int is_beg, int is_end)
{
  // Process new hits for a range of seeds.

  const int is_num = is_end - is_beg;

  // printf("CandCloner::ProcessSeedRange is_beg=%d, is_end=%d, is_num=%d\n", is_beg, is_end, is_num);

  //1) sort the candidates
  for (int is = is_beg; is < is_end; ++is)
  {
    auto& hitsForSeed = m_hits_to_add[is];
    auto& cands = mp_etabin_of_comb_candidates->m_candidates;

#ifdef DEBUG
    int th_start_seed = m_start_seed;

    std::cout << "dump seed n " << is << " with input candidates=" << hitsForSeed.size() << std::endl;
    for (int ih = 0; ih<hitsForSeed.size(); ih++)
    {
      std::cout << "trkIdx=" << hitsForSeed[ih].trkIdx << " hitIdx=" << hitsForSeed[ih].hitIdx << " chi2=" <<  hitsForSeed[ih].chi2 << std::endl;
      std::cout << "original pt=" << cands[th_start_seed+is][hitsForSeed[ih].trkIdx].pT() << " " 
                << "nTotalHits="  << cands[th_start_seed+is][hitsForSeed[ih].trkIdx].nTotalHits() << " " 
                << "nFoundHits="  << cands[th_start_seed+is][hitsForSeed[ih].trkIdx].nFoundHits() << " " 
                << "chi2="        << cands[th_start_seed+is][hitsForSeed[ih].trkIdx].chi2() << " " 
                << std::endl;
    }
#endif

    if ( ! hitsForSeed.empty())
    {
      int num_hits = std::min((int) hitsForSeed.size(), Config::maxCandsPerSeed);
      auto& cv = t_cands_for_next_lay[is - is_beg];

      for (int ih = 0; ih < num_hits; ih++)
      {
        MkFitter::IdxChi2List h2a = hitsForSeed.top();
        hitsForSeed.pop();
        auto& newcv = cv.emplace_start(cands[ m_start_seed + is ][ h2a.trkIdx ]);
        newcv.addHitIdx(h2a.hitIdx, 0);
        newcv.setChi2(h2a.chi2);
        cv.emplace_finish();
      }

      // Copy the best -2 cands back to the current list.
      if (num_hits < Config::maxCandsPerSeed)
      {
        std::vector<Track> &ov = cands[m_start_seed + is];
        int cur_m2 = 0;
        int max_m2 = ov.size();
        while (cur_m2 < max_m2 && ov[cur_m2].getLastHitIdx() != -2) ++cur_m2;
        while (cur_m2 < max_m2 && num_hits < Config::maxCandsPerSeed)
        {
          cv.maybe_push( ov[cur_m2++] );
          ++num_hits;
        }
      }

      cands[ m_start_seed + is ].swap(cv.rep());
      cv.clear();
    }
    // else
    // {
    //   // MT: make sure we have all cands with last hit idx == -2 at this point
    //
    //   for (auto &cand : cands[ m_start_seed + is ])
    //   {
    //     assert(cand.getLastHitIdx() == -2);
    //   }
    // }
  }
}

void CandCloner::DoWorkInSideThread(CandClonerWork_t work)
{
  int beg     = work.first;
  int the_end = work.second;

  // printf("CandCloner::DoWorkInSideThread working on beg=%d to the_end=%d\n", beg, the_end);

  while (beg != the_end)
  {
    int end = std::min(beg + s_max_seed_range, the_end);

    // printf("CandCloner::DoWorkInSideThread processing %4d -> %4d\n", beg, end);

    ProcessSeedRange(beg, end);

    beg = end;
  }
}
