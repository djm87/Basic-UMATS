           
      include 'modules_main.for'
      
      program main
      
      use miscellaneous_sub_main, only : tnorm,svdcmp,jacobi,solve_proc
     #    ,statev_inc
      use dimensions
      
      integer :: ibc(NTENS,10)
      real :: bc(NTENS,10),STATEV(NSTATV,10),DROT(3,3,10)
     #    ,DFGRD0(3,3,10),STRESS(NTENS,10),STRAN(NTENS,10),time(2,10)
     #    ,TEMP(10)
      character(len=6) :: str
      character(len=100) :: filename
      
      !apply bc from 0 state
      open(unit=111,file='bc.in')
      read(111,*)npro,NOELT
      read(111,*)tolerance
      
      do NOEL=1,NOELT
          write(str,'(I2)')NOEL
          lenpid=len_trim(str)
          filename='bc'//str(1:lenpid)//'.in'
          open(unit=111+NOEL,file=filename,status='unknown')
          filename='output'//str(1:lenpid)//'.out'
          open(unit=222+NOEL,file=filename,status='unknown')
      enddo
             
      do ipro=1,npro
          do NOEL=1,NOELT
              !define bc & time
              do i=1,3
                  DROT(i,i,NOEL)=1.0
                  DFGRD0(i,i,NOEL)=1.0
              enddo
              read(111+NOEL,*) (ibc(i,NOEL),i=1,NTENS)
              read(111+NOEL,*) (bc(i,NOEL),i=1,NTENS)
              read(111+NOEL,*) TSTART,TEND
              read(111+NOEL,*) time_pro,dtime
                        
              !solution for process ipro and element NOEL
              call solve_proc(NOEL,tolerance,ibc(:,NOEL),bc(:,NOEL),DROT
     #            ,DFGRD0,TSTART,TEND,time_pro,dtime,STATEV(:,NOEL)
     #            ,STRESS(:,NOEL),STRAN(:,NOEL),time(:,NOEL),TEMP(NOEL))
          enddo
      enddo
      
      close(111)
      do NOEL=1,NOELT
          close(111+NOEL)
          close(222+NOEL)
      enddo
      
      !apply bc from non 0 state (one inc)
c      call statev_inc
      
      end program main
      
      include 'UMAT.f'