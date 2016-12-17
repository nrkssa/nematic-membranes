	module module_energy_calculations
	Contains
!---------------------------------------------------------------------------------------------------------------------------
!                        SUBROUTINE TO CALCULATE THE ENERGY DIFFERENCE FOR A LINK FLIP
!---------------------------------------------------------------------------------------------------------------------------
         SUBROUTINE energy_linkflip(cver,tn1,tn2,energy)
         USE module_datastruct
         IMPLICIT NONE
         INTEGER,DIMENSION(4) :: cver
         INTEGER :: i,ip1,v,v1,tn1,tn2
         REAL(KIND=8) :: necuen,energy,spenergy
         energy=0.0
         flip_vertices: DO i=1,4,1
         necuen=0 ; spenergy=0 
         ip1=i+1 ; IF(i .Eq. 4) ip1=1
         v=cver(i)  ; v1=cver(ip1)                  									            ! Vertex involved are given by cver     
         Call Nematic_Interaction_Potential(spenergy,v,v1)       		          		    ! Nematic-Nematic interaction part      
		 Call Nematic_Curvature_Potential(necuen,v)
         energy=energy+(ver(v)%mcur**2)*kappa*ver(v)%totarea+spenergy+necuen                        !  Helfrich+Nem-nem+Nem-curvature       
         ENDDO flip_vertices
		 energy = energy-pr*(tri(tn1)%vol + tri(tn2)%vol)
         END SUBROUTINE energy_linkflip  

!------------------------------------------------------------------------------------------------------------------------
!                               SUBROUTINE TO CALCULATE THE ENERGY FOR VERTEX MOVES
!------------------------------------------------------------------------------------------------------------------------
         SUBROUTINE energy_vertexmove(vert,energy)
         USE module_datastruct
         IMPLICIT NONE
         INTEGER :: v,i,vert,v1,ip1,tn1 
         REAL(KIND=8) :: energy,necuen,spenergy
         energy=0.0

         mv_neigh: DO i=1,ver(vert)%nonei,1
         ip1=i+1 ; IF(i.EQ.ver(vert)%nonei) ip1=1
		 tn1 = ver(vert)%vneitr(i)
         necuen=0 ; spenergy=0.0
         v=ver(vert)%vneipt(i) ; v1=ver(vert)%vneipt(ip1) 
         Call Nematic_Interaction_Potential(spenergy,v,v1)                                    ! Nematic - Nematic interaction with all neigh spins   
	 	 Call Nematic_Curvature_Potential(necuen,v)
         energy=energy+(ver(v)%mcur**2)*kappa*ver(v)%totarea+spenergy+necuen-pr*tri(tn1)%vol          !  Helfrich favouring a positive spontaneous curv       
         ENDDO mv_neigh 
         necuen=0.0 ; Call Nematic_Curvature_Potential(necuen,vert)
         energy=energy+(ver(vert)%mcur**2)*kappa*ver(vert)%totarea+necuen                             !  Helfrich +Nematic-curvature for the displaced vertex 
         END SUBROUTINE energy_vertexmove

!------------------------------------------------------------------------------------------------------------------------
!                            SUBROUTINE TO CALCULATE THE ENERGY FOR MOVING THE NEMATIC VECTOR
!------------------------------------------------------------------------------------------------------------------------
          SUBROUTINE Nematic_Interaction_Potential(spenergy,v,v1)
          USE module_datastruct ; Use module_curvcalc; Use OMP_LIB
          IMPLICIT NONE
          INTEGER :: i,j,v
          INTEGER,OPTIONAL :: v1                                                              ! v1 is to skipped for overcounting in flips and moves
          REAL(KIND=8) :: tang,spenergy,nemint                                                ! The value of cos(\theta) and temproary storage      
          REAL(KIND=8),DIMENSION(3,3):: tapl1,tapl2                                           ! tangent planes 1 and 2                              
          REAL(KIND=8),DIMENSION(3,1):: geod,pg1,pg2,cp1,cp2,n1,n2                            ! Geodesic, projections 1,2, cross prod 1,2           
          REAL(KIND=8),DIMENSION(1,1):: c1,c2,s1,s2
            
          spenergy=0

          !$OMP Parallel Do Reduction(+:spenergy) Private(j,nemint,n1,n2,geod,tapl1,tapl2,pg1,pg2,c1,c2,s1,s2,cp1,cp2,tang)
          DO i=1,ver(v)%nonei,1
             j=ver(v)%vneipt(i) 
		     nemint=nempar(ver(v)%phase,ver(j)%phase)                                      ! The kind of nematic interaction between v and j   

          IF((Present(v1)) .And. (j.EQ.v1)) Then                                           ! Skipping the specified point 
		     spenergy = spenergy + 0.0
	      Else
              n1=ver(v)%vnor ; n2=ver(j)%vnor
              geod=ver(j)%vcoord-ver(v)%vcoord                                              ! Geodesic from vertex 1 to 2
              geod=geod/SQRT(SUM(geod**2))             
              tapl1=unitmat-MATMUL(n1,TRANSPOSE(n1))                                        ! Tangent plane at vertex v
              tapl2=unitmat-MATMUL(n2,TRANSPOSE(n2))                                        ! Tangent plane at vertex j
              pg1=MATMUL(tapl1,geod)                                                        ! Projection of geod on tang plane 
              pg1=pg1/SQRT(SUM(pg1**2))
              pg2=MATMUL(tapl2,geod)                                    
              pg2=pg2/SQRT(SUM(pg2**2))
              c1=MATMUL(TRANSPOSE(pg1),ver(v)%spgl)                                          ! cosine of angle of spin1 with the geodesic 
              c2=MATMUL(TRANSPOSE(pg2),ver(j)%spgl)                                          ! cosine of angle of spin2 with the geodesic 

              cp1(1,1)=n1(2,1)*pg1(3,1)-n1(3,1)*pg1(2,1)
              cp1(2,1)=n1(3,1)*pg1(1,1)-n1(1,1)*pg1(3,1) 
              cp1(3,1)=n1(1,1)*pg1(2,1)-n1(2,1)*pg1(1,1) 
              cp1=cp1/SQRT(SUM(cp1**2)) 
              s1=MATMUL(TRANSPOSE(cp1),ver(v)%spgl)                                         ! Sine of the angle between spin1 and geodesic  
              
              cp2(1,1)=n2(2,1)*pg2(3,1)-n2(3,1)*pg2(2,1)
              cp2(2,1)=n2(3,1)*pg2(1,1)-n2(1,1)*pg2(3,1) 
              cp2(3,1)=n2(1,1)*pg2(2,1)-n2(2,1)*pg2(1,1) 
              cp2=cp2/SQRT(SUM(cp2**2))
              s2=MATMUL(TRANSPOSE(cp2),ver(j)%spgl)                                         ! sine of angle between spin2 and the geodesic          
              tang=0.5*(3*(c1(1,1)*c2(1,1)+s1(1,1)*s2(1,1))**2-1)
              spenergy = spenergy -nemint*tang
	      	Endif
		ENDDO
	  	!$OMP End Parallel Do
        END SUBROUTINE Nematic_Interaction_Potential


!------------------------------------------------------------------------------------------------------------------------
!                       SUBROUTINE FOR COMPUTING THE NEMATIC CURVTURE INTERACTION AT A GIVEN VERTEX
!------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE Nematic_Curvature_Potential(nemenergy,v)
        USE module_datastruct 
        IMPLICIT NONE
        INTEGER :: v
		Real(Kind=8):: nemenergy,kpar,kper
		
        nemenergy=0.0; 	kpar = ver(v)%kpar; kper = ver(v)%kper
	
        If(kpar.Ne.0)Then                                                                              ! Do this only when this parameter is non zero                              
        nemenergy=nemenergy+kpar*((ver(v)%splo(1,1)**2*ver(v)%cur1+ver(v)%splo(2,1)**2*ver(v)%cur2)-ver(v)%cpar)**2 
        ! $\kappa_{\parallel}(nKn-c_0^{\parallel})^2$ (term along the orinetation)  
        Endif
        If(ver(v)%kper.Ne.0)Then 
        nemenergy=nemenergy+ver(v)%kper*((ver(v)%splo(2,1)**2*ver(v)%cur1+ver(v)%splo(1,1)**2*ver(v)%cur2)-ver(v)%cper)**2          
        ! Perpendicular interaction term (N is perpendicular to n)                  
        Endif
 		nemenergy=nemenergy*ver(v)%totarea                                                                                          
 		! All terms are integrated over the total area around the vertex  
        END SUBROUTINE Nematic_Curvature_Potential

!------------------------------------------------------------------------------------------------------------------------
!                       SUBROUTINE FOR COMPUTING THE NEMATIC CURVTURE INTERACTION AT A GIVEN VERTEX
!------------------------------------------------------------------------------------------------------------------------
        SUBROUTINE Helfrich_Potential(helfenergy,v)
        USE module_datastruct 
        IMPLICIT NONE
        INTEGER :: v
		Real(Kind=8)::helfenergy
        helfenergy = 0.0
		helfenergy = kappa*ver(v)%totarea*(ver(v)%mcur**2)
        END SUBROUTINE Helfrich_Potential

!------------------------------------------------------------------------------------------------------------------------
!                            SUBROUTINE TO CALCULATE THE ENERGY FOR MOVING THE NEMATIC VECTOR
!------------------------------------------------------------------------------------------------------------------------
        Function Nematic_Interaction_energy(v) Result(spenergy)
        IMPLICIT NONE
        INTEGER :: v
        REAL(KIND=8) :: spenergy                                                                ! The value of cos(\theta) and temproary storage      
        Call Nematic_Interaction_potential(spenergy,v)
        Return
        END function Nematic_Interaction_Energy

!------------------------------------------------------------------------------------------------------------------------
!                       SUBROUTINE FOR COMPUTING THE NEMATIC CURVTURE INTERACTION AT A GIVEN VERTEX
!------------------------------------------------------------------------------------------------------------------------
        Function Nematic_Curvature_Energy(v) Result(nemenergy)
        IMPLICIT NONE
        INTEGER :: v
		Real(Kind=8)::nemenergy
		Call Nematic_Curvature_Potential(nemenergy,v)
		Return
        END Function Nematic_Curvature_Energy

!------------------------------------------------------------------------------------------------------------------------
!                       SUBROUTINE FOR COMPUTING THE NEMATIC CURVTURE INTERACTION AT A GIVEN VERTEX
!------------------------------------------------------------------------------------------------------------------------
        Function Helfrich_Energy(v) Result(helfenergy)
        IMPLICIT NONE
        INTEGER :: v
		Real(Kind=8)::helfenergy
		Call Helfrich_Potential(helfenergy,v)
        END function Helfrich_Energy
	
	End Module module_energy_calculations
