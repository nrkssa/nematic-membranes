!===========================================================================================================================
!===========================================================================================================================
!				MODULE TO START FROM AN ARBITRARY CONFIGURATION
!===========================================================================================================================
!===========================================================================================================================

       MODULE module_rerun
        Integer,Allocatable,Dimension(:) :: nbm                        
        Integer,Allocatable,Dimension(:,:) :: ip
       contains 
       
!----------------------------------------------------------------------------------------------------------------------------
!                           Make a hexahedron and triangulate its surface
!----------------------------------------------------------------------------------------------------------------------------     
        SUBROUTINE make_hexahedron(N,b)
        USE module_datastruct
        Real(Kind=8) :: dl,theta,H,b
        Real(Kind=8),Allocatable,Dimension(:) :: x,y,z                     !ip --> a temp array to keep track of neigh   
        Integer::ipn,in,nt2,np,lz,ls
        Integer::ite,ll,N,nt,i
        nver=3*(N**2)+2                                                  !Number of vertices
        ntr=(nver-2)*2                                                   !Number of triangles
        tlink=nver+ntr-2                                                 ! Total number of links  
        ipn=0
        nt=3*N*(N+1)/2+1

        Allocate(x(nver))                                                !The size of the arrays involved are allocated 
        Allocate(y(nver))
        Allocate(z(nver))
        Allocate(ip(nver,9))
        Allocate(nbm(nver))  

        Allocate(ver(nver)) ; Allocate(tri(ntr));  Allocate(lin(-tlink:tlink))

         x=0; y=0; z=0; ip=0; nbm=0

         dl=N*b
         theta=asin(1.0/sqrt(3.0))
         H=dl*cos(theta)

               up_half_face1: do lz=0,N                                   !First face of the upper half plane
                  if(lz.lt.1)ll=0
                  if(lz.ge.1)ll=lz-1
                  do ls=0,ll
                  ipn=ipn+1
                  z(ipn)=H-b*lz*cos(theta)
                  x(ipn)=b*lz*sin(theta)/2.0
                  y(ipn)=-dl/2.0+ls*b+(N-lz)*b/2.0
                  enddo
               enddo up_half_face1
       
               up_half_face2: do lz=1,N                                   !Second face of the upper half plane
                  do ls=0,lz-1
                  ipn=ipn+1
                  z(ipn)=H-b*lz*cos(theta)
                  x(ipn)=b*(lz-ls)*sqrt(3.0)/2.0-b*lz*sin(theta)
                  y(ipn)=b*lz/2.0-ls*b/2.0
                  enddo
               enddo up_half_face2
           
               up_half_face3: do lz=1,N                                   ! Third face of the upper half plane
                     do ls=0,lz-1
                     ipn=ipn+1
                     z(ipn)=H-b*lz*cos(theta)
                     x(ipn)=b*ls*sqrt(3.0)/2.0-b*lz*sin(theta)
                     y(ipn)=-ls*b/2.0
                     enddo
               enddo up_half_face3

                   nbm(1)=3      
                   nbm(nt+1)=3                                           ! The number of neighbours for the top
                                                                         ! and bottom particle is 3
         unei_loop1: do i=1,3
                      ip(1,i)=(i-1)*N*(N+1)/2+2                          ! Neighbours of particle 1 is fixed in 
                                                                         ! every loop (2,17,32 for l=5)         
                       lz=1 ; ls=0
                       np=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls+2               ! The points nearest to the top point 
                       nbm(np)=6                                         ! All of these have 6 neigh each      

                       in=i+1                                
                       if(in.gt.3)in=1

                       ip(np,1)=(in-1)*N*(N+1)/2+2
                       ip(np,2)=1

                       in=i-1
                       if(in.lt.1)in=3

                       ip(np,3)=(in-1)*N*(N+1)/2+2
                       ip(np,4)=(in-1)*N*(N+1)/2+2+2
                       ip(np,5)=(i-1)*N*(N+1)/2+lz*(lz+1)/2+2
                       ip(np,6)=(i-1)*N*(N+1)/2+lz*(lz+1)/2+ls+1+2
        

           unei_loop2: do lz=2,N
             unei_loop3: do ls=0,lz-1
                         np=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls+2
                         nbm(np)=6

                  if(ls.ne.0.and.lz.ne.N)then
                       if(ls.ne.lz-1)then
                       ip(np,1)=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls+1+2
                       ip(np,2)=(i-1)*N*(N+1)/2+(lz-1)*(lz-2)/2+ls+2
                       else
                       ite=i*N*(N+1)/2+lz*(lz-1)/2+2
                       if(ite.gt.nt)ite=ite-nt+1
                       ip(np,1)=ite
                       ite=i*N*(N+1)/2+(lz-1)*(lz-2)/2+2
                       if(ite.gt.nt)ite=ite-nt+1
                       ip(np,2)=ite
                       endif

                       ip(np,3)=(i-1)*N*(N+1)/2+(lz-1)*(lz-2)/2+ls-1+2
                       ip(np,4)=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls-1+2
                       ip(np,5)=(i-1)*N*(N+1)/2+lz*(lz+1)/2+ls+2
                       ip(np,6)=(i-1)*N*(N+1)/2+lz*(lz+1)/2+ls+1+2

                  elseif(ls.eq.0.and.lz.ne.N)then
                      in=i-1
                      if(in.lt.1)in=3

                      ip(np,1)=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls+1+2
                      ip(np,2)=(i-1)*N*(N+1)/2+(lz-1)*(lz-2)/2+ls+2
                      ip(np,3)=(in-1)*N*(N+1)/2+lz*(lz-1)/2+lz-1+2
                      ip(np,4)=(in-1)*N*(N+1)/2+lz*(lz+1)/2+lz+2
                      ip(np,5)=(i-1)*N*(N+1)/2+lz*(lz+1)/2+ls+2
                      ip(np,6)=(i-1)*N*(N+1)/2+lz*(lz+1)/2+ls+1+2

                  elseif(ls.ne.0.and.lz.eq.N)then
                     if(ls.ne.lz-1)then
                     ip(np,1)=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls+1+2
                     ip(np,2)=(i-1)*N*(N+1)/2+(lz-1)*(lz-2)/2+ls+2
                     ip(np,6)=(i-1)*N*(N-1)/2+(N-1)*(N-2)/2+ls+2+nt
                     else
                     ite=i*N*(N+1)/2+lz*(lz-1)/2+2
                     if(ite.gt.nt)ite=ite-nt+1
                     ip(np,1)=ite
                     ite=i*N*(N+1)/2+(lz-1)*(lz-2)/2+2
                     if(ite.gt.nt)ite=ite-nt+1
                     ip(np,2)=ite
                     in=i+1
                     if(in.gt.3)in=1
                     ip(np,6)=(in-1)*N*(N-1)/2+(N-1)*(N-2)/2+2+nt
                     endif

                     ip(np,3)=(i-1)*N*(N+1)/2+(lz-1)*(lz-2)/2+ls-1+2
                     ip(np,4)=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls-1+2
                     ip(np,5)=(i-1)*N*(N-1)/2+(N-1)*(N-2)/2+ls-1+2+nt

                  elseif(ls.eq.0.and.lz.eq.N)then
                     in=i-1
                     if(in.lt.1)in=3
                     ip(np,1)=(i-1)*N*(N+1)/2+lz*(lz-1)/2+ls+1+2
                     ip(np,2)=(i-1)*N*(N+1)/2+(lz-1)*(lz-2)/2+ls+2
                     ip(np,3)=(in-1)*N*(N+1)/2+lz*(lz-1)/2+lz-1+2
                     ip(np,4)=(i-1)*N*(N-1)/2+(N-1)*(N-2)/2+2+nt
                     nbm(np)=4
                  endif
                enddo unei_loop3
             enddo  unei_loop2
           enddo   unei_loop1


                low_half_face1: do lz=0,N-1
                          if(lz.lt.1)ll=0
                          if(lz.ge.1)ll=lz-1
                           do ls=0,ll
                           ipn=ipn+1
                           z(ipn)=-H+b*lz*cos(theta)
                           x(ipn)=b*lz*sin(theta)/2.0
                           y(ipn)=-dl/2.0+ls*b+(N-lz)*b/2.0
                           enddo
                enddo low_half_face1

                low_half_face2: do lz=1,N-1
                        do ls=0,lz-1
                        ipn=ipn+1
                        z(ipn)=-H+b*lz*cos(theta)
                        x(ipn)=b*(lz-ls)*sqrt(3.0)/2.0-b*lz*sin(theta)
                        y(ipn)=b*lz/2.0-ls*b/2.0
                        enddo
                enddo low_half_face2

                low_half_face3: do lz=1,N-1
                     do ls=0,lz-1
                     ipn=ipn+1
                     z(ipn)=-H+b*lz*cos(theta)
                     x(ipn)=b*ls*sqrt(3.0)/2.0-b*lz*sin(theta)
                     y(ipn)=-ls*b/2.0
                     enddo
                enddo low_half_face3

                  nver=ipn 
                 
           lnei_loop1: do i=1,3
                       ip(nt+1,4-i)=(i-1)*N*(N-1)/2+2+nt
                       lz=1 ; ls=0
                       np=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+2+nt
                       nbm(np)=6

                       in=i+1
                       if(in.gt.3)in=1

                       ip(np,6)=(in-1)*N*(N-1)/2+2+nt
                       ip(np,5)=1+nt

                       in=i-1
                       if(in.lt.1)in=3

                       ip(np,4)=(in-1)*N*(N-1)/2+2+nt
                       ip(np,3)=(in-1)*N*(N-1)/2+4+nt
                       ip(np,2)=(i-1)*N*(N-1)/2+lz*(lz+1)/2+2+nt
                       ip(np,1)=(i-1)*N*(N-1)/2+lz*(lz+1)/2+ls+3+nt
                       nt2=3*N*(N-1)/2+1

            lnei_loop2: do lz=2,N-1
              lnei_loop3: do ls=0,lz-1
                      np=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+2+nt
                      nbm(np)=6

                   if(ls.ne.0.and.lz.ne.N-1)then
                      if(ls.ne.lz-1)then
                      ip(np,6)=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+3+nt
                      ip(np,5)=(i-1)*N*(N-1)/2+(lz-1)*(lz-2)/2+ls+2+nt
                      else
                      ite=i*N*(N-1)/2+lz*(lz-1)/2+2
                      if(ite.gt.nt2)ite=ite-nt2+1
                      ip(np,6)=ite+nt
                      ite=i*N*(N-1)/2+(lz-1)*(lz-2)/2+2
                      if(ite.gt.nt2)ite=ite-nt2+1
                      ip(np,5)=ite+nt
                      endif
                      ip(np,4)=(i-1)*N*(N-1)/2+(lz-1)*(lz-2)/2+ls+1+nt 
                      ip(np,3)=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+1+nt
                      ip(np,2)=(i-1)*N*(N-1)/2+lz*(lz+1)/2+ls+2+nt
                      ip(np,1)=(i-1)*N*(N-1)/2+lz*(lz+1)/2+ls+3+nt

                  elseif(ls.eq.0.and.lz.ne.N-1)then
                      in=i-1
                      if(in.lt.1)in=3
                      ip(np,6)=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+3+nt
                      ip(np,5)=(i-1)*N*(N-1)/2+(lz-1)*(lz-2)/2+ls+2+nt
                      ip(np,4)=(in-1)*N*(N-1)/2+lz*(lz-1)/2+lz+1+nt
                      ip(np,3)=(in-1)*N*(N-1)/2+lz*(lz+1)/2+lz+2+nt
                      ip(np,2)=(i-1)*N*(N-1)/2+lz*(lz+1)/2+ls+2+nt
                      ip(np,1)=(i-1)*N*(N-1)/2+lz*(lz+1)/2+ls+3+nt

                 elseif(ls.ne.0.and.lz.eq.N-1)then
                      if(ls.ne.lz-1)then
                      ip(np,6)=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+3+nt
                      ip(np,5)=(i-1)*N*(N-1)/2+(lz-1)*(lz-2)/2+ls+2+nt
                      else
                      ite=i*N*(N-1)/2+lz*(lz-1)/2+2
                      if(ite.gt.nt2)ite=ite-nt2+1
                      ip(np,6)=ite+nt
                      ite=i*N*(N-1)/2+(lz-1)*(lz-2)/2+2
                      if(ite.gt.nt2)ite=ite-nt2+1
                      ip(np,5)=ite+nt
                      endif
                      ip(np,4)=(i-1)*N*(N-1)/2+(lz-1)*(lz-2)/2+ls+1+nt
                      ip(np,3)=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+1+nt
                      ip(np,2)=(i-1)*N*(N+1)/2+N*(N-1)/2+ls+2
                      ip(np,1)=(i-1)*N*(N+1)/2+N*(N-1)/2+ls+3

                 elseif(ls.eq.0.and.lz.eq.N-1)then
                      in=i-1
                      if(in.lt.1)in=3
                      ip(np,6)=(i-1)*N*(N-1)/2+lz*(lz-1)/2+ls+3+nt
                      ip(np,5)=(i-1)*N*(N-1)/2+(lz-1)*(lz-2)/2+ls+2+nt
                      ip(np,4)=(in-1)*N*(N-1)/2+lz*(lz-1)/2+lz+1+nt
                      ip(np,3)=(in-1)*N*(N+1)/2+N*(N-1)/2+N+1
                      ip(np,2)=(i-1)*N*(N+1)/2+N*(N-1)/2+2
                      ip(np,1)=(i-1)*N*(N+1)/2+N*(N-1)/2+3
                 endif

                enddo lnei_loop3
             enddo lnei_loop2
           enddo lnei_loop1

          coord:DO i=1,nver                                              !Feed vertex coord and coord no to DS 
          ver(i)%vcoord=reshape((/x(i),y(i),z(i)/),shape(ver(i)%vcoord))
          ENDDO  coord 
          

          DEAllocate(x) ; DEAllocate(y) ; DEAllocate(z)

          Call order_triangles()
          
          Do i=1,ntr,1
           Call areacal(i)
          Enddo
          
          END SUBROUTINE make_hexahedron
!---------------------------------------------------------------------------------------------------------------------------
!                                         SUBROUTINE TO INDEX THE TRIANGLES AND LINKS
!---------------------------------------------------------------------------------------------------------------------------
         SUBROUTINE order_triangles()                                    !To order the vertices that make the triangles
          USE module_datastruct
          IMPLICIT NONE
          Integer:: kn,ikn,knp1,n2,nt,n1
          Integer::nlin,trno
          Integer,Allocatable,Dimension(:,:)::tl,trl
          Allocate(tl(nver,nver))
          Allocate(trl(nver,nver))
          tl=0 ; nt=0 ; nlin=0 ; trl=0


         over_allvert: DO kn=1,nver
           ver(kn)%nonei=nbm(kn)                                         !Store the number of neighbours to a vertex 
           ver(kn)%vneipt(1:nbm(kn))=ip(kn,1:nbm(kn))

          over_neigh: DO ikn=1,nbm(kn)                                   !To fix the coord number of the vertex
           knp1=ikn+1 ; IF(ikn.EQ.nbm(kn)) knp1=1
           n1=ip(kn,ikn) ; n2=ip(kn,knp1)


            IF(trl(kn,n1).EQ.0)THEN
            nt=nt+1
            tri(nt)%vert=reshape((/kn,n1,n2/),shape(tri(nt)%vert))
            trl(kn,n1)=nt ; trl(n1,n2)=nt; trl(n2,kn)=nt
            ver(kn)%vneitr(ikn)=nt
            call areacal(nt)

             IF(tl(kn,n1).EQ.0)THEN  
             nlin=nlin+1
             lin(nlin)%sep=reshape((/kn,n1/),shape(lin(nlin)%sep))
             lin(-nlin)%sep=reshape((/n1,kn/),shape(lin(-nlin)%sep))
             tl(kn,n1)=nlin ; tl(n1,kn)=-nlin
             ENDIF

             IF(tl(n1,n2).EQ.0)THEN
             nlin=nlin+1
             lin(nlin)%sep=reshape((/n1,n2/),shape(lin(nlin)%sep))
             lin(-nlin)%sep=reshape((/n2,n1/),shape(lin(-nlin)%sep))
             tl(n1,n2)=nlin ; tl(n2,n1)=-nlin 
             ENDIF

             IF(tl(n2,kn).EQ.0)THEN 
             nlin=nlin+1
             lin(nlin)%sep=reshape((/n2,kn/),shape(lin(nlin)%sep))
             lin(-nlin)%sep=reshape((/kn,n2/),shape(lin(-nlin)%sep))
             tl(n2,kn)=nlin ; tl(kn,n2)=-nlin 
             ENDIF
         
             lin(tl(kn,n1))%tr=nt ; lin(tl(n1,n2))%tr=nt 
             lin(tl(n2,kn))%tr=nt
             tri(nt)%li=(/tl(kn,n1),tl(n1,n2),tl(n2,kn)/) 

             ELSE
             trno=trl(kn,n1)
             lin(tl(kn,n1))%tr=trno
             lin(tl(n1,n2))%tr=trno
             lin(tl(n2,kn))%tr=trno
             ver(kn)%vneitr(ikn)=trno
             ENDIF 
               
             ENDDO over_neigh
           ENDDO over_allvert
           
           ntr=nt
           tlink=nlin
           DEAllocate(ip)
           DEAllocate(nbm)  
           DEAllocate(tl)

           END SUBROUTINE order_triangles         
!---------------------------------------------------------------------------------------------------------------------------
!                 Subroutine to calculate the area of the given face and update the total area linked to the vertices
!---------------------------------------------------------------------------------------------------------------------------

       SUBROUTINE areacal(tr)
       USE module_datastruct 
       IMPLICIT NONE
       Integer::i,j,k,tr
       REAL(KIND=8) ::ax,ay,az,area
       REAL(KIND=8),Dimension(3,1)::r1,r2,r3,r21,r31
       ax=0 ; ay=0 ; az=0 ; area=0

       i=tri(tr)%vert(1) ; j=tri(tr)%vert(2) ; k=tri(tr)%vert(3)
       r1=ver(i)%vcoord ; r2=ver(j)%vcoord ; r3=ver(k)%vcoord                                                  ! The coordinates of the 3 vertices      
       r21=r2-r1 ; r31=r3-r1                                                                                   ! Relative position vectors for area calc
       ver(i)%totarea=ver(i)%totarea-tri(tr)%ar/3.0                                                            ! Contrib of exisiting area to totalarea 
       ver(j)%totarea=ver(j)%totarea-tri(tr)%ar/3.0
       ver(k)%totarea=ver(k)%totarea-tri(tr)%ar/3.0

        ax=(r21(2,1)*r31(3,1)-r21(3,1)*r31(2,1))*0.5
        ay=(r21(3,1)*r31(1,1)-r21(1,1)*r31(3,1))*0.5
        az=(r21(1,1)*r31(2,1)-r21(2,1)*r31(1,1))*0.5
        area=SQRT(ax**2+ay**2+az**2)
        tri(tr)%ar=area                                       
        ax=ax/area ; ay=ay/area ;  az=az/area                                                                  ! Normalized area vector 

       tri(tr)%fnor=RESHAPE((/ax,ay,az/),SHAPE(tri(tr)%fnor)) 
       ver(i)%totarea=ver(i)%totarea+tri(tr)%ar/3.0                                                            ! Contrib of new area to the total area of vert
       ver(j)%totarea=ver(j)%totarea+tri(tr)%ar/3.0
       ver(k)%totarea=ver(k)%totarea+tri(tr)%ar/3.0

       tri(tr)%vol=(r1(1,1)*(r2(2,1)*r3(3,1)-r2(3,1)*r3(2,1))+ r1(2,1)*(r2(3,1)*r3(1,1)-r2(1,1)*r3(3,1))+ &
                    r1(3,1)*(r2(1,1)*r3(2,1)-r2(2,1)*r3(1,1)))/6.0


       END SUBROUTINE areacal     
!--------------------------------------------------------------------------------------------------------------------------
!                                     SUBROUTINE TO MAKE A RESTART FILE
!--------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE Write_Memb_Configuration(time)
        USE module_datastruct
        Integer:: i,time,Fin
        Character(100) :: nname,fname,tname

318    FORMAT(10(I4,1X))
320    FORMAT(F9.6,1X,F9.6,1X,F9.6,1X,I2)

        Write(nname,*),mynod
        Write(tname,*),time
        fname="./rundir-"//Trim(Adjustl(nname))//'/startdet-'//Trim(Adjustl(tname))//'.dat'
        Fin=11+mynod

	  OPEN(Fin,FILE=Trim(Adjustl(fname)))
          WRITE(Fin,*),nver,"vertexdet"   
          DO i=1,nver,1
          WRITE(Fin,*)ver(i)%vcoord
          WRITE(Fin,*)ver(i)%nonei
          WRITE(Fin,318)ver(i)%vneipt
          WRITE(Fin,318)ver(i)%vneitr
          ENDDO

          WRITE(Fin,*),ntr,"triangledet"
          DO i=1,ntr
          WRITE(Fin,*)tri(i)%vert
          WRITE(Fin,*)tri(i)%li
          ENDDO
   
          WRITE(Fin,*)tlink,"linkdet"
          DO i=-tlink,tlink,1
             IF(i.NE.0)THEN
             WRITE(FIn,*)lin(i)%sep  
             WRITE(Fin,*)lin(i)%tr  
             ENDIF
          ENDDO

          Write(Fin,*),nver,"nematicdet"      
          DO i=1,nver,1
          WRITE(Fin,320)ver(i)%splo,ver(i)%phase
          ENDDO

          Write(Fin,*),nver,"nematicdet-glo"      
          DO i=1,nver,1
          WRITE(Fin,320)ver(i)%spgl,ver(i)%phase
          ENDDO

          Write(Fin,*),nver,"c1  c2 mcur"      
          DO i=1,nver,1
          WRITE(Fin,*)ver(i)%cur1,ver(i)%cur2,ver(i)%mcur
          ENDDO

          Call MPI_Barrier(nallgrp,ier)	  
          CLOSE(Fin)

        END SUBROUTINE Write_Memb_Configuration

!---------------------------------------------------------------------------------------------------------------------------
!                              To start  from one compact startupfile
!---------------------------------------------------------------------------------------------------------------------------

       SUBROUTINE Read_Memb_Configuration(state)
       USE module_datastruct 
       Integer::i
       Character(100) :: state,nname,fname

		Write(nname,*),mynod
		fname='./startdet.in'
		OPEN(11,FILE=Trim(Adjustl(fname)),FORM='FORMATTED')
        READ(11,*),nver                                                ! Number of vertex
	  	Allocate(ver(nver))
          DO i=1,nver,1
          READ(11,*)ver(i)%vcoord                                        ! Position of coordinates 
          READ(11,*)ver(i)%nonei                                         ! Position of coordinates 
          READ(11,*)ver(i)%vneipt                                        ! Neighbouring vertex     
          READ(11,*)ver(i)%vneitr                                        ! Neighbouring triangles  
          ENDDO

          READ(11,*),ntr
		  Allocate(tri(ntr)) 
          DO i=1,ntr
          READ(11,*)tri(i)%vert
          READ(11,*)tri(i)%li
          ENDDO

          READ(11,*)tlink
		  Allocate(lin(-tlink:tlink))
          DO i=-tlink,tlink,1
             IF(i.NE.0)THEN
             READ(11,*)lin(i)%sep
             READ(11,*)lin(i)%tr
             ENDIF
          ENDDO

          If(Trim(state).Eq.'RESTART')Then                               ! Read the nematic data only if called for 
          READ(11,*)  
          DO i=1,nver,1
          READ(11,*)ver(i)%splo,ver(i)%phase
          ENDDO
          ElseIf(Trim(state).Eq.'MEMBRANE')Then                          ! Do not read the nematic data otherwise
          ver(1:nver)%splo(1,1)=0.0
          ver(1:nver)%splo(2,1)=0.0
          ver(1:nver)%splo(3,1)=0.0
          Endif

          CLOSE(11)

          DO i=1,ntr,1                                                   ! Calculate the area of the triangles 
          Call areacal(i)
          ENDDO
       
          END SUBROUTINE Read_Memb_Configuration

        END MODULE module_rerun

