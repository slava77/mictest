// HOST / AVX
// icc -O3 -mavx -I. -o mult66 mult66.cc  -vec-report=3
// ./mult66

// MIC
// icc -O3 -mmic -I. -o mult66-mic mult66.cc -vec-report=3
// scp mult66-mic root@mic0:
// ssh root@mic0 ./mult66-mic

#include "MatriplexSymNT.h"

#include "mplex-common.h"


typedef Matriplex<float, M, M> MPlexMM;
typedef MatriplexSym<float, M> MPlexSS;


int main(int arg, char *argv[])
{
  SMatrixSS *mul  = new_sth<SMatrixSS>(Ns);
  SMatrixSS *muz  = new_sth<SMatrixSS>(Ns);
  SMatrixMM *res  = new_sth<SMatrixMM>(Ns);

  MPlexSS mpl(N);
  MPlexSS mpz(N);
  MPlexMM mpres(N);
  

  init_mulmuz(mul, muz);

  // Multiplex arrays

  for (int i = 0; i < N; ++i)
  {
    mpl.Assign(i, mul[i].Array());
    mpz.Assign(i, muz[i].Array());
  }

  double t0, tmp, tsm;

  // ================================================================

  if (NN_MULT)
  {
    t0 = dtime();
    for (int i = 0; i < NN_MULT; ++i)
    {
      //#pragma omp simd collapse(8)
#pragma ivdep
      for (int m = 0; m < N; ++m)
      {
        res[m] = mul[m] * muz[m];
      }
    }

    tsm = dtime() - t0;
    std::cout << "SMatrix multiply time = " << tsm << " s\n";

    // ----------------------------------------------------------------

    t0 = dtime();

    for (int i = 0; i < NN_MULT; ++i)
    {
      Multiply(mpl, mpz, mpres);
    }


    tmp = dtime() - t0;
    std::cout << "Matriplex multiply time = " << tmp << " s\n";
    std::cout << "SMatrix / Matriplex = " << tsm/tmp << "";
    {
      double x = 0, y = 0;
      for (int j = 0; j < N; ++j)
      {
        x += res[j](1,2);
        y += mpres(1,2,j);
      }
      std::cout << "\t\t\tx = " << x << ", y = " << y << "\n";
    }
    std::cout << "\n";

    /*
      for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
      for (int k = 0; k < M; ++k)
      {
      if (res[i](j,k) != mpres.At(j, k, i))
      std::cout << i <<" "<< j <<" "<< k <<" "<< res[i](j,k) <<" "<< mpres.At(j, k, i) << "\n";
      }
    */
  }

  // ================================================================

  if (NN_INV)
  {
    t0 = dtime();

    for (int i = 0; i < NN_INV; ++i)
    {
      //#pragma omp simd collapse(8)
#pragma ivdep
      for (int m = 0; m < N; ++m)
      {
        mul[m].InvertFast();
        muz[m].InvertFast();
        //bool bl = mul[m].InvertChol();
        //bool bz = muz[m].InvertChol();

        //if ( ! bl || ! bz)   printf("Grr %d %d %d\n", m, bl, bz);
      }
    }
    tsm = dtime() - t0;
    std::cout << "SMatrix invert time = " << tsm << " s\n";

    // ----------------------------------------------------------------

    t0 = dtime();

    for (int i = 0; i < NN_INV; ++i)
    {
      SymInvertCramer(mpl);
      SymInvertCramer(mpz);
      // SymInvertChol(mpl);
      // SymInvertChol(mpz);
    }

    tmp = dtime() - t0;
    std::cout << "Matriplex invert time = " << tmp << " s\n";
    std::cout << "SMatrix / Matriplex = " << tsm/tmp << "";
    {
      double x = 0, y = 0;
      for (int j = 0; j < N; ++j)
      {
        x += mul[j](1,2) + muz[j](1,2);
        y += mpl(1,2,j)  + mpl(1,2,j);
      }
      std::cout << "\t\t\tx = " << x << ", y = " << y << "\n";
    }
    std::cout << "\n";

    // ================================================================

    if (COMPARE_INVERSES)
    {
      for (int i = 0; i < M; ++i) { for (int j = 0; j < M; ++j)
          printf("%8f ", mul[0](i,j)); printf("\n");
      } printf("\n");

      for (int i = 0; i < M; ++i) { for (int j = 0; j < M; ++j)
          printf("%8f ", mpl.At(i,j,0)); printf("\n");
      } printf("\n");
    }
  }

  return 0;
}