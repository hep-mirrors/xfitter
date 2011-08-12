      program tolhapdf
      real*8 buffer(8,0:160,0:160)
      real*8 q2valpdf(0:160),xvalpdf(0:160), alphas(0:160)
C------------------------------


      open (51,file='output/lhapdf.block.txt',status='old')
C First read Q2:
      do iq2=1,23
         read (51,'(7E12.4)')  (q2valpdf((iq2-1)*7+j),j=0,6)
      enddo
      do jx=1,23
         read (51,'(7E12.4)')  (xvalpdf((jx-1)*7+j),j=0,6)
      enddo
      do iq2=1,23
         read (51,'(7E12.4)')  (alphas((iq2-1)*7+j),j=0,6)
      enddo


      do i=0,160
         do j=0,160
            read (51,*) (buffer(k,j,i),k=1,8)
         enddo
      enddo
      close (51)

      open (51,file='lhapdf_tail.dat',status='unknown')

      do iq2=1,23
         write (51,'(7E12.4)') (q2valpdf((iq2-1)*7+j),j=0,6)
      enddo

      do jx=1,23
         write (51,'(7E12.4)') (xvalpdf((jx-1)*7+j),j=0,6)
      enddo

      do iq2=1,23
         write (51,'(7E12.4)') (alphas((iq2-1)*7+j),j=0,6)
      enddo


      do k=1,8
         do i=0,160
            do J=1,23
               write (51,'(7E12.4)')
     $              (buffer(k,(j-1)*7+ii-1,i),ii=1,7)
c               if (i.eq.0) then
c                  print *,(j-1)*7+ii
c               endif
            enddo            
         enddo
      enddo
      close (51)
      end
