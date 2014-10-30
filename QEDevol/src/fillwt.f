C     ------------------------------------------------------------------
C     Weight filling routines
C     ------------------------------------------------------------------

C     ==================================================================
      subroutine FilWTqcd(w,jset,idpij,idaij,nordqcd)
C     ==================================================================

C     Fill Pij unpolarised tables upto NNLO QCD

C     w           (in)   store previously partitioned by maketab
C     jset        (in)   set identifier
C     idpij       (out)  list of Pij table identifiers (global format)
C     idaij       (out)  list of Aij table identifiers (global format)
C     nordqcd     (in)   maximum QCD perturbative order (1 - LO, 2 - NLO, 3 - NNLO)

      implicit double precision (a-h,o-z)

      dimension w(*)
      dimension idpij(28,4)
      dimension idaij(5)

      external dqcAchi
      external dqcP111R,dqcP111S,dqcP111D
      external dqcP131A
      external dqcP231A
      external dqcP321A
      external dqcP331A,dqcP331R,dqcP331S,dqcP331D
      external dqcP112A,dqcP112B
      external dqcP122A,dqcP122B
      external dqcP132A
      external dqcP222A,dqcP222B
      external dqcP232A
      external dqcP322A
      external dqcP332A,dqcP332B
      external dqcP552B
      external dqcP113A,dqcP113B,dqcP113D
      external dqcP123A,dqcP123B,dqcP123D
      external dqcP133A
      external dqcP223A,dqcP223B,dqcP223D
      external dqcP233A
      external dqcP323A
      external dqcP333A,dqcP333B,dqcP333D
      external dqcP553A,dqcP553B,dqcP553D
      external dqcP563A
      external dqcP663A,dqcP663B,dqcP663D
      external dqcAQQ2A,dqcAQQ2B,dqcAQQ2D
      external dqcAGQ2A
      external dqcAGG2A,dqcAGG2B,dqcAGG2D
      external dqcAHQ2A
      external dqcAHG2A,dqcAHG2D

      write(6,'(/'' FilWTqcd fill user weight tables'')')

C     Loop over spline order
      do l = 1,3
        do k = 1,28
          idPij(k,l) = 1000*jset+200+28*(l-1)+k
        enddo
      enddo

      do i = 1,5
        idAij(i) = 1000*jset+100+i
      enddo

C     Pij LO QCD
      write(6,'('' FilWt Pij LO QCD'')')
      call MakeWRS(w,1000*jset+201,dqcP111R,dqcP111S,dqcAchi,0)
      call MakeWtD(w,1000*jset+201,dqcP111D,dqcAchi)
      call scalewt(w, 0.D0, 1000*jset+202)
      call MakeWtA(w,1000*jset+203,dqcP131A,dqcAchi)
      call scalewt(w,0.D0,1000*jset+204)
      call scalewt(w,0.D0,1000*jset+205)
      call CopyWgt(w,1000*jset+201,1000*jset+206,0)
      call MakeWtA(w,1000*jset+207,dqcP231A,dqcAchi)
      call scalewt(w,0.D0,1000*jset+208)
      call scalewt(w,0.D0,1000*jset+209)
      call MakeWtA(w,1000*jset+210,dqcP321A,dqcAchi)
      call MakeWtA(w,1000*jset+211,dqcP331A,dqcAchi)
      call MakeWRS(w,1000*jset+211,dqcP331R,dqcP331S,dqcAchi,0)
      call MakeWtD(w,1000*jset+211,dqcP331D,dqcAchi)
      call scalewt(w,0.D0,1000*jset+212)
      call scalewt(w,0.D0,1000*jset+213)
      call scalewt(w,0.D0,1000*jset+214)
      call scalewt(w,0.D0,1000*jset+215)
      call scalewt(w,0.D0,1000*jset+216)
      call CopyWgt(w,1000*jset+201,1000*jset+217,0)
      call scalewt(w, 0.D0, 1000*jset+218)
      call scalewt(w, 0.D0, 1000*jset+219)
      call CopyWgt(w,1000*jset+201,1000*jset+220,0)
      call CopyWgt(w,1000*jset+201,1000*jset+221,0)
      call CopyWgt(w,1000*jset+201,1000*jset+222,0)
      call CopyWgt(w,1000*jset+201,1000*jset+223,0)
      call CopyWgt(w,1000*jset+201,1000*jset+224,0)
      call CopyWgt(w,1000*jset+201,1000*jset+225,0)
      call CopyWgt(w,1000*jset+201,1000*jset+226,0)
      call CopyWgt(w,1000*jset+201,1000*jset+227,0)
      call CopyWgt(w,1000*jset+201,1000*jset+228,0)

C     Pij NLO QCD
      write(6,'('' FilWt Pij NLO QCD'')')
      if (nordqcd.lt.2) then
        do k = 1,28
          call scalewt(w,0.D0,1000*jset+228+k)
        enddo
      else
        call MakeWtA(w,1000*jset+229,dqcP112A,dqcAchi)
        call MakeWtB(w,1000*jset+229,dqcP112B,dqcAchi,0)
        call MakeWtA(w,1000*jset+230,dqcP122A,dqcAchi)
        call MakeWtB(w,1000*jset+230,dqcP122B,dqcAchi,0)
        call MakeWtA(w,1000*jset+231,dqcP132A,dqcAchi)
        call scalewt(w,0.D0,1000*jset+232)
        call scalewt(w,0.D0,1000*jset+233)
        call MakeWtA(w,1000*jset+234,dqcP222A,dqcAchi)
        call MakeWtB(w,1000*jset+234,dqcP222B,dqcAchi,0)
        call MakeWtA(w,1000*jset+235,dqcP232A,dqcAchi)
        call scalewt(w,0.D0,1000*jset+236)
        call scalewt(w,0.D0,1000*jset+237)
        call MakeWtA(w,1000*jset+238,dqcP322A,dqcAchi)
        call MakeWtA(w,1000*jset+239,dqcP332A,dqcAchi)
        call MakeWtB(w,1000*jset+239,dqcP332B,dqcAchi,0)
        call scalewt(w,0.D0,1000*jset+240)
        call scalewt(w,0.D0,1000*jset+241)
        call scalewt(w,0.D0,1000*jset+242)
        call scalewt(w,0.D0,1000*jset+243)
        call scalewt(w,0.D0,1000*jset+244)
        call MakeWtB(w,1000*jset+245,dqcP552B,dqcAchi,0)
        call scalewt(w,0.D0,1000*jset+246)
        call scalewt(w,0.D0,1000*jset+247)
        call CopyWgt(w,1000*jset+245,1000*jset+248,0)
        call CopyWgt(w,1000*jset+229,1000*jset+249,0)
        call CopyWgt(w,1000*jset+229,1000*jset+250,0)
        call CopyWgt(w,1000*jset+229,1000*jset+251,0)
        call CopyWgt(w,1000*jset+229,1000*jset+252,0)
        call CopyWgt(w,1000*jset+245,1000*jset+253,0)
        call CopyWgt(w,1000*jset+245,1000*jset+254,0)
        call CopyWgt(w,1000*jset+245,1000*jset+255,0)
        call CopyWgt(w,1000*jset+245,1000*jset+256,0)
      endif

C     Pij NNLO QCD
      write(6,'('' FilWt Pij NNLO QCD'')')
      if (nordqcd.lt.3) then
        do k = 1,28
          call scalewt(w,0.D0,1000*jset+256+k)
        enddo
      else
        call MakeWtA(w,1000*jset+257,dqcP113A,dqcAchi)
        call MakeWtB(w,1000*jset+257,dqcP113B,dqcAchi,1)
        call MakeWtD(w,1000*jset+257,dqcP113D,dqcAchi)
        call MakeWtA(w,1000*jset+258,dqcP123A,dqcAchi)
        call MakeWtA(w,1000*jset+259,dqcP133A,dqcAchi)
        call scalewt(w,0.D0,1000*jset+260)
        call scalewt(w,0.D0,1000*jset+261)
        call MakeWtA(w,1000*jset+262,dqcP223A,dqcAchi)
        call MakeWtB(w,1000*jset+262,dqcP223B,dqcAchi,1)
        call MakeWtD(w,1000*jset+262,dqcP223D,dqcAchi)
        call MakeWtA(w,1000*jset+263,dqcP233A,dqcAchi)
        call scalewt(w,0.D0,1000*jset+264)
        call scalewt(w,0.D0,1000*jset+265)
        call MakeWtA(w,1000*jset+266,dqcP323A,dqcAchi)
        call MakeWtA(w,1000*jset+267,dqcP333A,dqcAchi)
        call MakeWtB(w,1000*jset+267,dqcP333B,dqcAchi,1)
        call MakeWtD(w,1000*jset+267,dqcP333D,dqcAchi)
        call scalewt(w,0.D0,1000*jset+268)
        call scalewt(w,0.D0,1000*jset+269)
        call scalewt(w,0.D0,1000*jset+270)
        call scalewt(w,0.D0,1000*jset+271)
        call scalewt(w,0.D0,1000*jset+272)
        call MakeWtA(w,1000*jset+273,dqcP553A,dqcAchi)
        call MakeWtB(w,1000*jset+273,dqcP553B,dqcAchi,1)
        call MakeWtD(w,1000*jset+273,dqcP553D,dqcAchi)
        call MakeWtA(w,1000*jset+274,dqcP563A,dqcAchi)
        call scalewt(w,0.D0,1000*jset+275)
        call MakeWtA(w,1000*jset+276,dqcP663A,dqcAchi)
        call MakeWtB(w,1000*jset+276,dqcP663B,dqcAchi,1)
        call MakeWtD(w,1000*jset+276,dqcP663D,dqcAchi)
        call CopyWgt(w,1000*jset+257,1000*jset+277,0)
        call CopyWgt(w,1000*jset+257,1000*jset+278,0)
        call CopyWgt(w,1000*jset+257,1000*jset+279,0)
        call CopyWgt(w,1000*jset+257,1000*jset+280,0)
        call CopyWgt(w,1000*jset+273,1000*jset+281,0)
        call CopyWgt(w,1000*jset+273,1000*jset+282,0)
        call CopyWgt(w,1000*jset+273,1000*jset+283,0)
        call CopyWgt(w,1000*jset+273,1000*jset+284,0)
      endif

C     Aij NNLO QCD
      write(6,'('' FilWt Aij NNLO QCD'')')
      if (nordqcd.lt.3) then
        do k = 1,5
          call scalewt(w,0.D0,1000*jset+100+k)
        enddo
      else
        call MakeWtA(w,1000*jset+101,dqcAQQ2A,dqcAchi)
        call MakeWtB(w,1000*jset+101,dqcAQQ2B,dqcAchi,1)
        call MakeWtD(w,1000*jset+101,dqcAQQ2D,dqcAchi)
        call MakeWtA(w,1000*jset+102,dqcAGQ2A,dqcAchi)
        call MakeWtA(w,1000*jset+103,dqcAGG2A,dqcAchi)
        call MakeWtB(w,1000*jset+103,dqcAGG2B,dqcAchi,1)
        call MakeWtD(w,1000*jset+103,dqcAGG2D,dqcAchi)
        call MakeWtA(w,1000*jset+104,dqcAHQ2A,dqcAchi)
        call MakeWtA(w,1000*jset+105,dqcAHG2A,dqcAchi)
        call MakeWtD(w,1000*jset+105,dqcAHG2D,dqcAchi)
      endif

      return
      end

C     ==================================================================
      subroutine FilWTqed(w,jset,idpij,nordqed)
C     ==================================================================

C     Fill Pij unpolarised tables at LO QED

C     w           (in)   store previously partitioned by maketab
C     jset        (in)   set identifier
C     idpij       (out)  list of Pij table identifiers (global format)
C     nordqed     (in)   maximum QED perturbative order (0 - no QED, 1 - LO QED)

      implicit double precision (a-h,o-z)

      dimension w(*)
      dimension idpij(28,4)

      external dqcAchi
      external dqcP114R,dqcP114S,dqcP114D
      external dqcP124R,dqcP124S,dqcP124D
      external dqcP144A
      external dqcP244A
      external dqcP414A
      external dqcP424A
      external dqcP444D
      external dqcP554R,dqcP554S,dqcP554D
      external dqcP564R,dqcP564S,dqcP564D
      external dqcP774R,dqcP774S,dqcP774D
      external dqcP884R,dqcP884S,dqcP884D

      write(6,'(/'' FilWTqed fill user weight tables'')')

C     Loop over spline order
      do k = 1,28
        idPij(k,4) = 1000*jset+200+k
      enddo

C     Pij LO QED
      write(6,'('' FilWt Pij LO QED'')') 
      if (nordqed.lt.1) then
        do k = 1,28
          call scalewt(w,0.D0,1000*jset+200+k)
        enddo
      else
        call MakeWRS(w,1000*jset+201,dqcP114R,dqcP114S,dqcAchi,0)
        call MakeWtD(w,1000*jset+201,dqcP114D,dqcAchi)
        call MakeWRS(w,1000*jset+202,dqcP124R,dqcP124S,dqcAchi,0)
        call MakeWtD(w,1000*jset+202,dqcP124D,dqcAchi)
        call scalewt(w,0.D0,1000*jset+203)
        call MakeWtA(w,1000*jset+204,dqcP144A,dqcAchi)
        call CopyWgt(w,1000*jset+202,1000*jset+205,0)
        call CopyWgt(w,1000*jset+201,1000*jset+206,0)
        call scalewt(w,0.D0,1000*jset+207)
        call MakeWtA(w,1000*jset+208,dqcP244A,dqcAchi)
        call scalewt(w,0.D0,1000*jset+209)
        call scalewt(w,0.D0,1000*jset+210)
        call scalewt(w,0.D0,1000*jset+211)
        call scalewt(w,0.D0,1000*jset+212)
        call MakeWtA(w,1000*jset+213,dqcP414A,dqcAchi)
        call MakeWtA(w,1000*jset+214,dqcP424A,dqcAchi)
        call scalewt(w,0.D0,1000*jset+215)
        call MakeWtD(w,1000*jset+216,dqcP444D,dqcAchi)
        call MakeWRS(w,1000*jset+217,dqcP554R,dqcP554S,dqcAchi,0)
        call MakeWtD(w,1000*jset+217,dqcP554D,dqcAchi)
        call MakeWRS(w,1000*jset+218,dqcP564R,dqcP564S,dqcAchi,0)
        call MakeWtD(w,1000*jset+218,dqcP564D,dqcAchi)
        call CopyWgt(w,1000*jset+218,1000*jset+219,0)
        call CopyWgt(w,1000*jset+217,1000*jset+220,0)
        call MakeWRS(w,1000*jset+221,dqcP774R,dqcP774S,dqcAchi,0)
        call MakeWtD(w,1000*jset+221,dqcP774D,dqcAchi)
        call MakeWRS(w,1000*jset+222,dqcP884R,dqcP884S,dqcAchi,0)
        call MakeWtD(w,1000*jset+222,dqcP884D,dqcAchi)
        call CopyWgt(w,1000*jset+221,1000*jset+223,0)
        call CopyWgt(w,1000*jset+222,1000*jset+224,0)
        call CopyWgt(w,1000*jset+221,1000*jset+225,0)
        call CopyWgt(w,1000*jset+222,1000*jset+226,0)
        call CopyWgt(w,1000*jset+221,1000*jset+227,0)
        call CopyWgt(w,1000*jset+222,1000*jset+228,0)
      endif

      return
      end
