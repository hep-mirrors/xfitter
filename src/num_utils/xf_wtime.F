C----------------------------------------------------------------------
C> @file xf_wtime.F
C> @brief Portable wall-clock timer.
C>
C> Returns wall-clock time as a double precision value (seconds).
C> When compiled with OpenMP (-fopenmp), uses omp_get_wtime.
C> Otherwise falls back to cpu_time, which reports CPU time rather than
C> wall-clock time but keeps timing prints functional.
C----------------------------------------------------------------------
      double precision function xf_wtime()
      implicit none
#ifdef _OPENMP
      double precision omp_get_wtime
      external omp_get_wtime
      xf_wtime = omp_get_wtime()
#else
      real t
      call cpu_time(t)
      xf_wtime = dble(t)
#endif
      return
      end
