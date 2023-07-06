#include <iostream>
#include <vector>
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/shm.h>
#include <tuple>
#include <algorithm>
#include <hf_errlog.h>

extern "C" {
  void fcall_(double* FSHERE, double* amin, double* x0best, const double (*FCN)(double*), const double (*FUTIL)(double*));
  void set_x_values_(const double* Xinput);
  double getumat_(int const& i, int const& j);
  void scaleumat_(int  const& i, int const& j, double const& scale);
  void compute_hessian_(const int& npar,
                        const int& NCPU,
                        const double* x0,
                        double* x0best,
                        const double& dstep,
                        const double& amin,
                        double& amin4,
                        double* vhvecOut,
                        double& cormax,
                        double& cormin,
                        double* G11Out,
                        const double (*FCN)(double*),
                        const double (*FUTIL)(double*));
  // has side effects on umat
  void step_method_one_(const int& npar,
                        const int& NCPU,
                        const double* x0,
                        double* x0best,
                        const double& aimsag,
                        const double& amin,
                        double& amin4,
                        const int& ncyc,
                        const double (*FCN)(double*),
                        const double (*FUTIL)(double*)
                       );
  // has side effects on umat
  void step_method_two_(const int& npar,
                        const int& NCPU,
                        const double* x0,
                        double* x0best,
                        const double& aimsag,
                        const double& amin,   //global
                        double& amin4,  //to be updated
                        const int& ncyc,
                        const double (*FCN)(double*),
                        const double (*FUTIL)(double*)
                       );
}


inline double umat(int i, int j)
{
  return getumat_(i + 1, j + 1);
}

inline void scaleUmat(int i, int j, double scale)
{
  scaleumat_(i + 1, j + 1, scale);
}

void compute_hessian_(const int& npar,
                      const int& NCPU,
                      const double* x0,
                      double* x0best,
                      const double& dstep,
                      const double& amin,
                      double& amin4,
                      double* vhvecOut,
                      double& cormax,
                      double& cormin,
                      double* G11Out,
                      const double (*FCN)(double*),
                      const double (*FUTIL)(double*))
{
  const bool debug = false;
  
  const double cmult = dstep / std::sqrt(2.0);
  int ntot = npar * (npar + 1) / 2;

  int shmid = shmget(IPC_PRIVATE, ntot * sizeof(double), IPC_CREAT | 0666);
  if (shmid < 0) {
    hf_errlog(202301501,"F: Shared memory allocation failed.");
  }
  double* vhvec = (double*)shmat(shmid, NULL, 0);

  int shmid2 = shmget(IPC_PRIVATE, ntot * sizeof(double), IPC_CREAT | 0666);
  if (shmid2 < 0) {
    hf_errlog(202301502,"F: Shared memory allocation failed.");
  }
  double* cmax = (double*)shmat(shmid2, NULL, 0);

  int shmid3 = shmget(IPC_PRIVATE, ntot * sizeof(double), IPC_CREAT | 0666);
  if (shmid3 < 0) {
    hf_errlog(202301503,"F: Shared memory allocation failed.");
  }
  double* cmin = (double*)shmat(shmid3, NULL, 0);

  int shmid4 = shmget(IPC_PRIVATE, npar * sizeof(double), IPC_CREAT | 0666);
  if (shmid4 < 0) {
    hf_errlog(202301504,"F: Shared memory allocation failed.");
  }
  double* G11 = (double*)shmat(shmid4, NULL, 0);


  // we also want to save accidently found better minimum
  int shmid5 = shmget(IPC_PRIVATE, (npar + 1) * NCPU * sizeof(double), IPC_CREAT | 0666);
  if (shmid5 < 0) {
    hf_errlog(202301505,"F: Shared memory allocation failed.");
  }
  // npar+1 to stor amin first
  double* x0bestCPU = (double*)shmat(shmid5, NULL, 0);

  double x[npar];

  // diagonal part, 2 calls per parameter. Loop over separately, as in iterate.
  int chunkSize = npar / NCPU;
  int reminder  = npar % NCPU;
  int startIndex = 0;
  int endIndex = 0;
  for (int cpu = 0; cpu < std::min(NCPU, npar); cpu++) {

    int startIndex = endIndex;
    endIndex   = startIndex + chunkSize;
    if (cpu < reminder) {
      endIndex += 1;
    }

    pid_t pid;
    pid = fork();

    if (pid == 0) {
      for (int ipar = startIndex; ipar < endIndex; ipar += 1) {
        // Take steps in opposite directions
        for (int j = 0; j < npar; j++) {
          x[j] = x0[j] + umat(j, ipar) * dstep;
        }
        double FS1;
        set_x_values_(x);
        fcall_(&FS1, &amin4, x0best, FCN, FUTIL);

        for (int j = 0; j < npar; j++) {
          x[j] = x0[j] - umat(j, ipar) * dstep;
        }
        double FS2;
        set_x_values_(x);
        fcall_(&FS2, &amin4, x0best, FCN, FUTIL);

        double SAG = 0.5 * (FS1 + FS2 - 2.0 * amin);
        double G22 = 2.0 * SAG / (dstep * dstep);   // Estimate of second derivative
        G11[ipar] = (FS1 - FS2) / (2.0 * dstep);       // Estimate of first derivative
        int NDEX = (ipar + 1) * (ipar + 2) / 2 - 1;
        vhvec[NDEX] = G22;
      }
      // here we can store x0
      x0bestCPU[cpu * (npar + 1) + 0] = amin4;
      for (int ipar = 0; ipar < npar; ipar += 1)
        x0bestCPU[cpu * (npar + 1) + 1 + ipar] = x0best[ipar];

      exit(0);
    } else if (pid < 0) {
      hf_errlog(202301520,"F: Fork failed.");
    }
  }

  int status;
  while (wait(&status) > 0);

  // store output for diagonal part
  for (int ipar = 0; ipar < npar; ipar += 1) {
    int NDEX = (ipar + 1) * (ipar + 2) / 2 - 1;
    vhvecOut[NDEX] = vhvec[NDEX];
    G11Out[ipar]   = G11[ipar];
  }

  // (potentially) update best minimum
  for (int cpu = 0; cpu < std::min(NCPU, npar); cpu += 1) {
    if (x0bestCPU[cpu * (npar + 1) + 0] < amin4) {
      if (debug) std::cout << " BETTER AMIN " << amin4 << " " << amin4 - x0bestCPU[cpu * (npar + 1) + 0] << std::endl;
      amin4 = x0bestCPU[cpu * (npar + 1) + 0];
      for (int ipar = 0; ipar < npar; ipar += 1) {
        if (debug) std::cout << " BETTER " << ipar << " " << x0best[ipar] << " " <<  x0bestCPU[cpu * (npar + 1) + 1 + ipar] << std::endl;
        x0best[ipar] = x0bestCPU[cpu * (npar + 1) + 1 + ipar];
      }
    }
  }


  // off-diagonal part (4 calls per pair,  (a,b) (-a,b) (a,-b) (-a,-b) )
  std::vector<int> idI(ntot);
  std::vector<int> idJ(ntot);
  std::vector<int> idN;
  for (int i = 2; i <= npar; i += 1)
    for (int j = 1; j < i; j += 1) {
      int idx = (i * (i - 1)) / 2 + j;
      idN.push_back(idx - 1);
      idI[idx - 1] = i - 1;
      idJ[idx - 1] = j - 1;
    };

  //  std::cout << ntot << " " <<  idN.size() << std::endl;
  ntot = idN.size();

  chunkSize = ntot / NCPU;
  reminder  = ntot % NCPU;
  startIndex = 0;
  endIndex = 0;

  for (int cpu = 0; cpu < std::min(NCPU, ntot); cpu++) {

    int startIndex = endIndex;
    endIndex   = startIndex + chunkSize;
    if (cpu < reminder) {
      endIndex += 1;
    }

    pid_t pid;
    pid = fork();

    if (pid == 0) {
      for (int idx = startIndex; idx < endIndex; idx += 1) {

        int NDEX = idN[idx];
        int I = idI[NDEX];
        int J = idJ[NDEX];

        //  std::cout << I << " " << J << " " << NDEX  << "\n";

        for (int k = 0; k < npar; k++) {
          x[k] = x0[k] + cmult * (umat(k, I) + umat(k, J));
        }
        double FS1;
        set_x_values_(x);
        fcall_(&FS1, &amin4, x0best, FCN, FUTIL);

        for (int k = 0; k < npar; k++) {
          x[k] = x0[k] - cmult * (umat(k, I) + umat(k, J));
        }
        double FS2;
        set_x_values_(x);
        fcall_(&FS2, &amin4, x0best, FCN, FUTIL);

        for (int k = 0; k < npar; k++) {
          x[k] = x0[k] + cmult * (umat(k, I) - umat(k, J));
        }
        double FS3;
        set_x_values_(x);
        fcall_(&FS3, &amin4, x0best, FCN, FUTIL);

        for (int k = 0; k < npar; k++) {
          x[k] = x0[k] - cmult * (umat(k, I) - umat(k, J));
        }
        double FS4;
        set_x_values_(x);
        fcall_(&FS4, &amin4, x0best, FCN, FUTIL);

        double slow = (FS1 + FS2 - FS3 - FS4) / (4.0 * cmult * cmult);
        //  std::cout << "NNN " << NDEX+1 << " " <<  " " << I+1 << " " << J+1 << " " << FS1-amin << " " << FS2-amin << " " << FS3-amin << " " << FS4-amin << std::endl;
        cmax[NDEX] = std::max({FS1 - amin4, FS2 - amin4, FS3 - amin4, FS4 - amin4});
        cmin[NDEX] = std::min({FS1 - amin4, FS2 - amin4, FS3 - amin4, FS4 - amin4});
        vhvec[NDEX] = slow;
      }
      // here we can store x0
      x0bestCPU[cpu * (npar + 1) + 0] = amin4;
      for (int ipar = 0; ipar < npar; ipar += 1)
        x0bestCPU[cpu * (npar + 1) + 1 + ipar] = x0best[ipar];

      exit(0);

    } else if (pid < 0) {
      hf_errlog(202301521,"F: Fork failed.");
    }
  }

  while (wait(&status) > 0);

  // (potentially) update best minimum
  for (int cpu = 0; cpu < std::min(NCPU, ntot); cpu += 1) {
    if (x0bestCPU[cpu * (npar + 1) + 0] < amin4) {
      if (debug) std::cout << " BETTER AMIN " << amin4 << " " << amin4 - x0bestCPU[cpu * (npar + 1) + 0] << std::endl;
      amin4 = x0bestCPU[cpu * (npar + 1) + 0];
      for (int ipar = 0; ipar < npar; ipar += 1) {
        if (debug) std::cout << " BETTER " << ipar << " " << x0best[ipar] << " " <<  x0bestCPU[cpu * (npar + 1) + 1 + ipar] << std::endl;
        x0best[ipar] = x0bestCPU[cpu * (npar + 1) + 1 + ipar];
      }
    }
  }

  cormax = cmax[0];
  cormin = cmin[0];
  for (int i = 0; i < ntot; i += 1) {
    int NDEX = idN[i];
    vhvecOut[NDEX] = vhvec[NDEX];
    if (cmax[i] > cormax)
      cormax = cmax[i];
    if (cmin[i] < cormin)
      cormin = cmin[i];
  }
  shmdt(vhvec);
  shmctl(shmid, IPC_RMID, NULL);
  shmdt(cmax);
  shmctl(shmid2, IPC_RMID, NULL);
  shmdt(cmin);
  shmctl(shmid3, IPC_RMID, NULL);
  shmdt(G11);
  shmctl(shmid4, IPC_RMID, NULL);
  shmdt(x0bestCPU);
  shmctl(shmid5, IPC_RMID, NULL);
}

// Other parts of the code


#include <math.h>

void step_method_one_(const int& npar,
                      const int& NCPU,
                      const double* x0,
                      double* x0best,
                      const double& aimsag,
                      const double& amin,
                      double& amin4,
                      const int& ncyc,
                      const double (*FCN)(double*),
                      const double (*FUTIL)(double*)
                     )
{
  const bool debug = false;
  
  double tmp, tmpst1, tmpst2, fs1, fs2, rat1, rat2;

  double x[npar];

  int shmid = shmget(IPC_PRIVATE, npar * sizeof(double), IPC_CREAT | 0666);
  if (shmid < 0) {
    hf_errlog(202301506,"F: Shared memory allocation failed.");
  }
  double* scaleShared = (double*)shmat(shmid, NULL, 0);

  // we also want to save accidently found better minimum
  int shmid2 = shmget(IPC_PRIVATE, (npar + 1) * NCPU * sizeof(double), IPC_CREAT | 0666);
  if (shmid2 < 0) {
    hf_errlog(202301507,"F: Shared memory allocation failed.");
  }
  // npar+1 to stor amin first
  double* x0bestCPU = (double*)shmat(shmid2, NULL, 0);

  int chunkSize = npar / NCPU;
  int reminder  = npar % NCPU;
  int startIndex = 0;
  int endIndex = 0;
  for (int cpu = 0; cpu < std::min(NCPU, npar); cpu++) {

    int startIndex = endIndex;
    endIndex   = startIndex + chunkSize;
    if (cpu < reminder) {
      endIndex += 1;
    }

    pid_t pid;
    pid = fork();

    if (pid == 0) {
      for (int ipar = startIndex; ipar < endIndex; ipar += 1) {
        double scale = 1.;
        for (int isize = 0; isize < ncyc; isize++) {
          for (int j = 0; j < npar; j++) {
            x[j] = x0[j] + umat(j, ipar) * sqrt(aimsag) * scale;
          }
          set_x_values_(x);
          fcall_(&fs1, &amin4, x0best, FCN, FUTIL);

          for (int j = 0; j < npar; j++) {
            x[j] = x0[j] - umat(j, ipar) * sqrt(aimsag) * scale;
          }
          set_x_values_(x);
          fcall_(&fs2, &amin4, x0best, FCN, FUTIL);

          rat1 = (fmax(fs1, fs2) - amin) / aimsag;
          rat2 = (fmin(fs1, fs2) - amin) / aimsag;

          if (rat2 <= 0.0) {
            tmp = 2.0;
            if (isize > 2) {
              if (tmpst1 > 1.0) tmp = fmin(tmp, tmpst1);
              if (tmpst2 > 1.0) tmp = fmin(tmp, tmpst2);
            }
          } else {
            tmp = sqrt(2.0 / (rat1 + rat2));
          }

          if (tmp > 1.5) {
            tmp = 2.0 - 0.75 / tmp;
          } else if (tmp < 0.7) {
            tmp = 0.5 + (2.0 / 7.0) * tmp;
          }

          if (isize > 2) {
            if (std::max({fabs(tmp - 1.0),
                          fabs(tmpst1 - 1.0),
                          fabs(tmpst2 - 1.0)}) < 0.1) {
              break;
            }
          }

          if (isize > 1) tmpst2 = tmpst1;
          tmpst1 = tmp;
          // Update umat
          scale *= tmp;
        }
        scaleShared[ipar] = scale;
      }
      // here we can store x0
      x0bestCPU[cpu * (npar + 1) + 0] = amin4;
      for (int ipar = 0; ipar < npar; ipar += 1)
        x0bestCPU[cpu * (npar + 1) + 1 + ipar] = x0best[ipar];
      exit(0);
    } else  if (pid < 0) {
      hf_errlog(202301522,"F: Fork failed.");
    }
  }
  int status;
  while (wait(&status) > 0);

  // scale umat
  for (int ipar = 0; ipar < npar; ipar += 1)
    for (int j = 0; j < npar; j++)
      scaleUmat(j, ipar, scaleShared[ipar]);
  // also update minimum, if better is found
  for (int cpu = 0; cpu < std::min(NCPU, npar); cpu += 1) {
    if (x0bestCPU[cpu * (npar + 1) + 0] < amin4) {
      if (debug) std::cout << " BETTER AMIN " << amin4 << " " << amin4 - x0bestCPU[cpu * (npar + 1) + 0] << std::endl;
      amin4 = x0bestCPU[cpu * (npar + 1) + 0];
      for (int ipar = 0; ipar < npar; ipar += 1) {
        if (debug) std::cout << " BETTER " << ipar << " " << x0best[ipar] << " " <<  x0bestCPU[cpu * (npar + 1) + 1 + ipar] << std::endl;
        x0best[ipar] = x0bestCPU[cpu * (npar + 1) + 1 + ipar];
      }
    }
  }
  shmdt(scaleShared);
  shmctl(shmid, IPC_RMID, NULL);
  shmdt(x0bestCPU);
  shmctl(shmid2, IPC_RMID, NULL);

}

void step_method_two_(const int& npar,
                      const int& NCPU,
                      const double* x0,
                      double* x0best,
                      const double& aimsag,
                      const double& amin,   //global
                      double& amin4,  //to be updated
                      const int& ncyc,
                      const double (*FCN)(double*),
                      const double (*FUTIL)(double*)
                     )
{
  bool const debug = false;
  double tmvec[2];
  double x[npar];

  int shmid = shmget(IPC_PRIVATE, npar * sizeof(double), IPC_CREAT | 0666);
  if (shmid < 0) {
    hf_errlog(202301508,"F: Shared memory allocation failed.");
  }
  double* scaleShared = (double*)shmat(shmid, NULL, 0);

  // we also want to save accidently found better minimum
  int shmid2 = shmget(IPC_PRIVATE, (npar + 1) * NCPU * sizeof(double), IPC_CREAT | 0666);
  if (shmid2 < 0) {
    hf_errlog(202301509,"F: Shared memory allocation failed.");
  }
  // npar+1 to stor amin first
  double* x0bestCPU = (double*)shmat(shmid2, NULL, 0);

  int chunkSize = npar / NCPU;
  int reminder  = npar % NCPU;
  int startIndex = 0;
  int endIndex = 0;
  for (int cpu = 0; cpu < std::min(NCPU, npar); cpu++) {

    int startIndex = endIndex;
    endIndex   = startIndex + chunkSize;
    if (cpu < reminder) {
      endIndex += 1;
    }

    pid_t pid;
    pid = fork();

    if (pid == 0) {
      for (int ipar = startIndex; ipar < endIndex; ipar += 1) {

        int igood = 0;
        double tmult;
        double xmult;
        for (int idir = 1; idir < 3; idir++) {
          tmult = sqrt(aimsag);
          if (idir == 2) {
            tmult = -tmult;
          }
          int ntry = 50;
          double dbesth = -9.99e99;
          double dbestl = 9.99e99;
          double tbesth = 1.0;
          double tbestl = 1.0;
          double tbest;
          double dbest;
          for (int itry = 0; itry < ntry; itry++) {
            for (int j = 0; j < npar; j++) {
              x[j] = x0[j] + umat(j, ipar) * tmult;
            }
            if (debug) std::cout << "itry " << itry << std::endl;

            double FS1;
            set_x_values_(x);
            fcall_(&FS1, &amin4, x0best, FCN, FUTIL);
            double del = FS1 - amin;


            if (debug) std::cout << "del " << del  << " fs1 " << FS1 << " amin " << amin << " tmult " << tmult << std::endl;

            double const accu = 1e-3;
            if (fabs(del - aimsag) < accu * aimsag) {
              igood += idir;
              std::cout << "all good" << aimsag << " " << std::endl;
              break;
            }

            if (itry == 0) {
              tbest = tmult;
              dbest = del;
            }

            if (del >= aimsag) {
              if (fabs(del - aimsag) < fabs(dbesth - aimsag)) {
                tbesth = tmult;
                dbesth = del;
              }
            }

            if (del <= aimsag) {
              if (fabs(del - aimsag) < fabs(dbestl - aimsag)) {
                tbestl = tmult;
                dbestl = del;
              }
            }

            if (debug) std::cout << "LH: tbesth " << tbesth << " dbesth " << dbesth << " tbestl " << tbestl <<  " dbetsl " << dbestl <<
                                   std::endl;

            if (del <= 0.0) {
              break;
            }

            if (itry == 0) {
              xmult = sqrt(aimsag / del);
              double ascale = 0.3;
              xmult = 1.0 + ascale * atan((xmult - 1.0) / ascale);
              goto label_850;
            }

            if (debug) std::cout << "AFTER LABEL tbesth " << tbesth << " dbesth " << dbesth << " tbestl " << tbestl <<  " dbetsl " << dbestl <<
                                   std::endl;

            if (fabs(del - aimsag) < fabs(dbest - aimsag)) {
              tbest = tmult;
              dbest = del;
            }

            if (((xmult > 1.0) && (del >= aimsag)) || ((xmult < 1.0) && (del <= aimsag))) {
              tmult = tbest;
              del = dbest;
              xmult = sqrt(xmult);
            }

            if (((xmult > 1.0) && (del >= aimsag)) || ((xmult < 1.0) && (del <= aimsag))) {
              xmult = 1.0 / xmult;
            }

label_850:
            tmult *= xmult;
            if (debug) std::cout << "at label 850" << xmult << " " << tmult << std::endl;

          }

          if ((igood != 3) && (idir > 1)) {
            tmult = tbest;

            for (int j = 0; j < npar; j++) {
              x[j] = x0[j] + umat(j, ipar) * tmult;
            }

            double FS1;
            set_x_values_(x);
            fcall_(&FS1, &amin4, x0best, FCN, FUTIL);

            printf("ITERATE: poor convergence warning -- Delta Chisqr=%.5e ipar=%d idir=%d\n", FS1 - amin, ipar, idir);
          }

          tmvec[idir - 1] = tmult;
        }

        if (igood == 3) {
          double tmp = 1.0 / (tmvec[0] * tmvec[0]) + 1.0 / (tmvec[1] * tmvec[1]) + 1.0 / (tmvec[0] * tmvec[1]);
          tmult = 1.0 / sqrt(tmp);
        } else if (igood == 1) {
          tmult = fabs(tmvec[0]);
        } else if (igood == 2) {
          tmult = fabs(tmvec[1]);
        } else {
          tmult = 0.001;
        }

        scaleShared[ipar] =  tmult / sqrt(aimsag);
        if (debug) std::cout << "TMUL " << ipar << " " << tmult
                               << " " << tmvec[0] - 1 << " " << tmvec[1] - 1 <<  std::endl;
      }
      // here we can store x0
      x0bestCPU[cpu * (npar + 1) + 0] = amin4;
      for (int ipar = 0; ipar < npar; ipar += 1)
        x0bestCPU[cpu * (npar + 1) + 1 + ipar] = x0best[ipar];
      exit(0);
    } else  if (pid < 0) {
      hf_errlog(202301523,"F: Fork failed.");
    }
  }
  int status;
  while (wait(&status) > 0);

  // scale umat
  for (int ipar = 0; ipar < npar; ipar += 1)
    for (int j = 0; j < npar; j++)
      scaleUmat(j, ipar, scaleShared[ipar]);
  // also update minimum, if better is found
  for (int cpu = 0; cpu < std::min(NCPU, npar); cpu += 1) {
    if (x0bestCPU[cpu * (npar + 1) + 0] < amin4) {
      if (debug) std::cout << " BETTER AMIN " << amin4 << " " << amin4 - x0bestCPU[cpu * (npar + 1) + 0] << std::endl;
      amin4 = x0bestCPU[cpu * (npar + 1) + 0];
      for (int ipar = 0; ipar < npar; ipar += 1) {
        if (debug) std::cout << " BETTER " << ipar << " " << x0best[ipar] << " " <<  x0bestCPU[cpu * (npar + 1) + 1 + ipar] << std::endl;
        x0best[ipar] = x0bestCPU[cpu * (npar + 1) + 1 + ipar];
      }
    }
  }
  shmdt(scaleShared);
  shmctl(shmid, IPC_RMID, NULL);
  shmdt(x0bestCPU);
  shmctl(shmid2, IPC_RMID, NULL);
}
