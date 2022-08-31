MODULE saw_def
  INTEGER, PARAMETER :: Width=12, Wh=(Width+1)/2
  INTEGER, PARAMETER :: NumThreads=4
  INTEGER(KIND=8), PARAMETER ::   HashMask=4**Wh-1, One=1, Two=2, Three=3
  INTEGER, PARAMETER :: Empty=0, Lower=1, Upper=2

  INTEGER(KIND=8) :: Motzkin_0(0:Width+5,-1:Width+7), Motzkin_1(0:Width+5,-1:Width+7)

  INTEGER(KIND=8) :: NumLeftSignatures(0:Width),  NumRightSignatures(0:Width)   ! Number of signatures of given height
  INTEGER(KIND=8) :: BaseLeftSignatures(0:Width),  BaseRightSignatures(0:Width) ! Base for storing signatures in single array
  INTEGER(KIND=8) :: Base(0:Width)                                              ! Base for hash function
  INTEGER(KIND=8), ALLOCATABLE :: LeftSignatures(:), RightSignatures(:)         ! Signatures
  INTEGER(KIND=8), ALLOCATABLE :: LeftHash(:), RightHash(:)                     ! Hash function for signatures coded as look-up table
  INTEGER(KIND=4), ALLOCATABLE :: TaskThread(:)                                 ! Allocation of signatures to threads coded as look-up table

  INTEGER(KIND=8), ALLOCATABLE :: WalkCount(:)    ! The number af crossing walks
  INTEGER(KIND=8) :: Prime

  INTEGER :: Prime15(0:29), Prime30(0:29), Prime62(0:29), Prime_Res, First_Prime, Last_Prime
  REAL(KIND=8) :: CurTime,CPUTime,WallTime

  DATA Prime15/19,49,51,55,61,75,81,115,121,135,147,157,159,165,181,189, &
       195,199,205,207,231,235,237,261,265,271,277,289,301,325/
  DATA Prime30/35,41,83,101,105,107,135,153,161,173,203,257,263,297, &
       321,347,357,383,405,425,437,443,453,495,513,515,537,587,611,627/
  DATA Prime62/57,87,117,143,153,167,171,195,203,273,287,317,443,483,495,&
       575,581,603,633,663,765,773,777,791,813,831,923,981,993,1001/

END MODULE saw_def

PROGRAM SAW

  !     Enumeration of honeycomb lattice self-avoiding walks crossing a rhombus patch.
  !
  !     Program last changed by Iwan Jensen 16-12-2021
  !


  USE saw_def
  USE OMP_LIB
  IMPLICIT NONE

  READ(*,*) First_Prime,Last_Prime
  CALL INITIALIZE
  DO Prime_Res=First_Prime,Last_Prime

     WallTime=OMP_GET_WTIME() 
     CALL CPU_TIME(CPUTime)

     CALL FIRSTCOLUMN
     CALL ADDCOLUMNS

     CurTime=OMP_GET_WTIME() 
     OPEN(7,FILE='hexsaw_cross_rhombus_omp_time.log',ACCESS='APPEND')
     WRITE(7,'(A,I3,A,I3,A,2F12.4)') "Width: ",Width,"  Threads: ",NumThreads,"  Wall Time: ",CurTime-WallTime
     WallTime=CurTime
     CALL CPU_TIME(CurTime)
     WRITE(7,'(A,I3,A,I3,A,2F12.4)') "Width: ",Width,"  Threads: ",NumThreads,"  & 
                                CPU Time (total and per thread):  ",CurTime-CPUTime,(CurTime-CPUTime)/NumThreads
     CPUTime=CurTime
     CLOSE(7)

  END DO

END PROGRAM SAW


SUBROUTINE INITIALIZE

  USE saw_def
  USE OMP_LIB
  IMPLICIT NONE
  INTEGER(KIND=8) :: K

  CALL OMP_SET_NUM_THREADS(NumThreads)


  CALL MotzkinNumbers()  ! Calculate number of Motzkin paths

  CALL ConstructHashTable(Wh)

  K=Motzkin_0(Width+2,0)-Motzkin_0(Width+1,0)
  ALLOCATE(WalkCount(K))            

  ALLOCATE(TaskThread(0:4**(Wh+1)))    ! Thread tables for parallel updating


END SUBROUTINE INITIALIZE


SUBROUTINE FIRSTCOLUMN

  !
  !	Start WITH top signature S0 = (0,0....,0,0,0): BP(S0) = 1 
  !

  USE saw_def
  IMPLICIT NONE
  INTEGER(KIND=8):: Sig,Sl,Sr,Ps

  Prime=2**30-Prime30(Prime_Res)
  Prime=Two**62-Prime62(Prime_Res)

  WalkCount=0

  Sig=ISHFT(One,2*Width+1)  ! Free end at top most edge
  Sl=IAND(Sig,HashMask)
  Sr=ISHFT(Sig,-2*Wh)
  Ps=LeftHash(Sl)+RightHash(Sr)
  WalkCount(Ps)=1

END SUBROUTINE FIRSTCOLUMN




SUBROUTINE ADDCOLUMNS

  !
  !	Put in a cell of the lattice:     
  !
  !
  !       |                                       |
  !     \ |       / TOP                   \       | /
  !      \|      /                         \      |/
  !       \     /                           \     /
  !       |\   /                             \   /|
  !       | \ /                               \ / |
  !       |  |                                 |  |
  !       |  | MIDDLE                          |  |
  !       |  |            ------------>        |  |
  !       |  |                                 |  |
  !       |  |                                 |  |
  !       | / \                               / \ |
  !       |/   \                             /   \|
  !       /     \                           /     \
  !      /|      \                         /      |\
  !     / |       \ BOTTOM                /       | \
  !       |                                       |
  !      TMB                                     TMB                        
  !      
  ! 
  !    Schematic illustation of how the transfer matrix boundary (TMB) is moved
  !    in order to insert a new 'cell'. 
  !    The TM updating is determined by the states of the TOP and BOTTOM edges intersected
  !    by the TMB prior to the move. Bonds can be inserted along the BOTTOM, MIDDLE and TOP
  !    edges as long as only valid partial SAW configurations are formed. 
  !    
  !
  !

  USE saw_def
  USE OMP_LIB
  IMPLICIT NONE
  INTEGER, PARAMETER ::  Wlt=Wh-1, Wlb=Wlt-1
  INTEGER :: I,J,K,M,KL,KU,Col,Row,Iam
  INTEGER(KIND=8) :: Sig,St,Sl,Sr,Ps,Pt,S1,S2
  REAL(KIND=8) :: PTime,CTime
  CHARACTER(LEN=5) :: CW
  CHARACTER(LEN=30) :: SFN


  !$OMP PARALLEL PRIVATE(Iam,Sig,St,Sl,Sr,Ps,Pt,S1,S2,Row,I,J,K,M,KL,KU,Ptime,Ctime)

  PTime=OMP_GET_WTIME() 
  Iam = omp_get_thread_num()

  DO Row=0,Width ! Build up patch column-by-column
     IF (Iam==0) THEN
        CALL ConstructSignatures(Wlt)
        CALL ConstructThreadTable(Wlt,Wlt)
     END IF
     !$OMP BARRIER

     DO KL=Width-1,Wlt,-1      ! Loop over the top half of a column
        KU=KL+1

        DO I=0,MIN(Wlt+1,Width+1-Wlt)                                     ! Loop over height of Left signature
           DO J=1,NumLeftSignatures(I)                                    ! Loop over Left signatures of height I
              Sl=LeftSignatures(J+BaseLeftSignatures(I))                  ! Left signature
              IF (Iam==TaskThread(Sl)) THEN                               ! Check to see which thread does the work      
                 DO K=1,NumRightSignatures(I)                             ! Loop over Right signatures of height I (independently for each thread)
                    Sr=RightSignatures(K+BaseRightSignatures(I))          ! Right signature
                    Sig=Sl+ISHFT(Sr,2*Wlt)                                ! Complete source signature by concatenation
                    S1=IAND(Sig,HashMask)                                 ! Left half of source signature for hashing
                    S2=ISHFT(Sig,-2*Wh)                                   ! Right half of source signature for hashing
                    Ps=LeftHash(S1)+RightHash(S2)                         ! Hash value where count is stored
                    
                    SELECT CASE (IAND(ISHFT(Sig,-2*KL),Three+4*Three))
                       
                    CASE (Empty+4*Empty)
                       
                       ! Bottom edge is in state Empty and Top edge is in state Empty

                       ! Insert partial new loop
                       St=Sig+ISHFT(One,2*KL)+ISHFT(One,2*KU+1)   ! Target signature with new arc
                       S1=IAND(St,HashMask)                       ! Left half of target signature for hashing
                       S2=ISHFT(St,-2*Wh)                         ! Right half of target signature for hashing
                       Pt=LeftHash(S1)+RightHash(S2)              ! Storage position of target count
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime)   ! Add it up

                    CASE (Lower+4*Empty)
                       
                       ! Bottom edge is in state Lower and Top edge is in state Empty
                       ! Note that (E,L) and (L,E) has the same count after update by symmetry

                       ! Continue Lower arc end along bottom edge or continue Lower arc end along top edge

                       St=Sig-ISHFT(One,2*KL)+ISHFT(One,2*KU)  ! Target signature with Lower arc end at Top edge
                       S1=IAND(St,HashMask)                        
                       S2=ISHFT(St,-2*Wh)                         
                       Pt=LeftHash(S1)+RightHash(S2)               
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime) 
                       WalkCount(Ps)=WalkCount(Pt)



                    CASE (Upper+4*Empty)
                       
                       ! Bottom edge is in state Upper and Top edge is in state Empty
                       ! Note that (E,U) and (U,E) has the same count after update by symmetry

                       ! Continue Upper arc end along Bottom edge or continue Upper arc end along Top edge

                       St=Sig-ISHFT(One,2*KL+1)+ISHFT(One,2*KU+1)   ! Target signature with Upper arc end at Top edge
                       S1=IAND(St,HashMask)           
                       S2=ISHFT(St,-2*Wh)
                       Pt=LeftHash(S1)+RightHash(S2)
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime)
                       WalkCount(Ps)=WalkCount(Pt)

                    CASE (Lower+4*Lower)

                       ! Bottom edge is in state Lower and Top edge is in state Lower

                       ! Join Lower arc ends and relabel matching Upper arc end as Lower
                       St=Sig-5*ISHFT(One,2*KL)
                       CALL Relabel_Upper(St,KL)
                       S1=IAND(St,HashMask)
                       S2=ISHFT(St,-2*Wh)
                       Pt=LeftHash(S1)+RightHash(S2)
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime)



                    CASE (Upper+4*Lower)
                       
                       ! Bottom edge is in state Upper and Top edge is in state Lower

                       ! Join the two arc ends  
                       St=Sig-6*ISHFT(One,2*KL)
                       S1=IAND(St,HashMask)
                       S2=ISHFT(St,-2*Wh)
                       Pt=LeftHash(S1)+RightHash(S2)
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime)

                    CASE (Upper+4*Upper)
                       
                       ! Bottom edge is in state Upper and Top edge is in state Upper

                       ! Join Upper arc ends and relabel matching Lower arc end as Upper
                       St=Sig-10*ISHFT(One,2*KL)
                       CALL Relabel_Lower(St,KL)
                       S1=IAND(St,HashMask)
                       S2=ISHFT(St,-2*Wh)
                       Pt=LeftHash(S1)+RightHash(S2)
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime)

                    END SELECT
                 END DO
              END IF
           END DO
        END DO
     END DO


     !$OMP BARRIER

     IF (Iam==0) THEN
        OPEN(7,FILE='hexsaw_cross_rhombus_omp.log',ACCESS='APPEND')
        CTime=OMP_GET_WTIME() 
        WRITE(7,'(A,I3,A,2F12.4)') "Calc -- Top: ",Row,"  Wall Time: ",CTime-Ptime
        PTime=Ctime
        CLOSE(7)
        CALL ClearSignatures

        CALL ConstructSignatures(Wlb)
     END IF

     !$OMP BARRIER

     KL=Wlt-1 
     KU=KL+1
     IF (Iam==0) CALL ConstructThreadTable(Wlb,KL)
     !$OMP BARRIER

     DO I=0,MIN(Wlb+1,Width+1-Wlb)                                 ! Loop over height of Left signature
        DO J=1,NumLeftSignatures(I)                                ! Loop over Left signatures of height I
           Sl=LeftSignatures(J+BaseLeftSignatures(I))              ! Left signature
           IF (Iam==TaskThread(Sl)) THEN                           ! Check to see which thread does the work  
              DO K=1,NumRightSignatures(I)                         ! Loop over Right signatures of height I
                 Sr=RightSignatures(K+BaseRightSignatures(I))      ! Right signature
                 Sig=Sl+ISHFT(Sr,2*Wlb)                            ! Complete source signature by concatenation
                 S1=IAND(Sig,HashMask)                             ! Left half of source signature for hashing
                 S2=ISHFT(Sig,-2*Wh)                               ! Right half of source signature for hashing
                 Ps=LeftHash(S1)+RightHash(S2)                     ! Hash value where count is stored
                    
                    SELECT CASE (IAND(ISHFT(Sig,-2*KL),Three+4*Three))
                       
                    CASE (Empty+4*Empty)
                       
                       ! Bottom edge is in state Empty and Top edge is in state Empty

                       ! Insert partial new loop
                       St=Sig+ISHFT(One,2*KL)+ISHFT(One,2*KU+1)   ! Target signature with new arc
                       S1=IAND(St,HashMask)                       ! Left half of target signature for hashing
                       S2=ISHFT(St,-2*Wh)                         ! Right half of target signature for hashing
                       Pt=LeftHash(S1)+RightHash(S2)              ! Storage position of target count
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime)   ! Add it up

                    CASE (Lower+4*Empty)
                       
                       ! Bottom edge is in state Lower and Top edge is in state Empty
                       ! Note that (E,L) and (L,E) has the same count after update by symmetry

                       ! Continue Lower arc end along bottom edge or continue Lower arc end along top edge

                       St=Sig-ISHFT(One,2*KL)+ISHFT(One,2*KU)  ! Target signature with Lower arc end at Top edge
                       S1=IAND(St,HashMask)                        
                       S2=ISHFT(St,-2*Wh)                         
                       Pt=LeftHash(S1)+RightHash(S2)               
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime) 
                       WalkCount(Ps)=WalkCount(Pt)



                    CASE (Upper+4*Empty)
                       
                       ! Bottom edge is in state Upper and Top edge is in state Empty
                       ! Note that (E,U) and (U,E) has the same count after update by symmetry

                       ! Continue Upper arc end along Bottom edge or continue Upper arc end along Top edge

                       St=Sig-ISHFT(One,2*KL+1)+ISHFT(One,2*KU+1)   ! Target signature with Upper arc end at Top edge
                       S1=IAND(St,HashMask)           
                       S2=ISHFT(St,-2*Wh)
                       Pt=LeftHash(S1)+RightHash(S2)
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime)
                       WalkCount(Ps)=WalkCount(Pt)

                    CASE (Lower+4*Lower)

                       ! Bottom edge is in state Lower and Top edge is in state Lower

                       ! Join Lower arc ends and relabel matching Upper arc end as Lower
                       St=Sig-5*ISHFT(One,2*KL)
                       CALL Relabel_Upper(St,KL)
                       S1=IAND(St,HashMask)
                       S2=ISHFT(St,-2*Wh)
                       Pt=LeftHash(S1)+RightHash(S2)
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime)



                    CASE (Upper+4*Lower)
                       
                       ! Bottom edge is in state Upper and Top edge is in state Lower

                       ! Join the two arc ends  
                       St=Sig-6*ISHFT(One,2*KL)
                       S1=IAND(St,HashMask)
                       S2=ISHFT(St,-2*Wh)
                       Pt=LeftHash(S1)+RightHash(S2)
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime)

                    CASE (Upper+4*Upper)
                       
                       ! Bottom edge is in state Upper and Top edge is in state Upper

                       ! Join Upper arc ends and relabel matching Lower arc end as Upper
                       St=Sig-10*ISHFT(One,2*KL)
                       CALL Relabel_Lower(St,KL)
                       S1=IAND(St,HashMask)
                       S2=ISHFT(St,-2*Wh)
                       Pt=LeftHash(S1)+RightHash(S2)
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime)

                    END SELECT
              END DO
           END IF
        END DO
     END DO

     !$OMP BARRIER


     IF (Iam==0) THEN
        OPEN(7,FILE='hexsaw_cross_rhombus_omp.log',ACCESS='APPEND')
        CTime=OMP_GET_WTIME() 
        WRITE(7,'(A,I3,A,2F12.4)') "Calc -- Mid: ",Row,"  Wall Time: ",CTime-Ptime
        PTime=Ctime
        CLOSE(7)
        CALL ClearSignatures

        CALL ConstructSignatures(Wlt)
     END IF

     !$OMP BARRIER

     DO KL=Wlt-2,0,-1     ! Loop over the bottom half of a column
        KU=KL+1
        IF (Iam==0) CALL ConstructThreadTable(Wlt,KL)
        !$OMP BARRIER

        DO I=0,MIN(Wlt+1,Width+1-Wlt)                                 ! Loop over height of Left signature
           DO J=1,NumLeftSignatures(I)                                ! Loop over Left signatures of height I
              Sl=LeftSignatures(J+BaseLeftSignatures(I))              ! Left signature
              IF (Iam==TaskThread(Sl)) THEN                           ! Check to see which thread does the work  
                 DO K=1,NumRightSignatures(I)                         ! Loop over Right signatures of height I
                    Sr=RightSignatures(K+BaseRightSignatures(I))      ! Right signature
                    Sig=Sl+ISHFT(Sr,2*Wlt)                            ! Complete source signature by concatenation
                    S1=IAND(Sig,HashMask)                             ! Left half of source signature for hashing
                    S2=ISHFT(Sig,-2*Wh)                               ! Right half of source signature for hashing
                    Ps=LeftHash(S1)+RightHash(S2)                     ! Hash value where count is stored
                    
                    SELECT CASE (IAND(ISHFT(Sig,-2*KL),Three+4*Three))
                       
                    CASE (Empty+4*Empty)
                       
                       ! Bottom edge is in state Empty and Top edge is in state Empty

                       ! Insert partial new loop
                       St=Sig+ISHFT(One,2*KL)+ISHFT(One,2*KU+1)   ! Target signature with new arc
                       S1=IAND(St,HashMask)                       ! Left half of target signature for hashing
                       S2=ISHFT(St,-2*Wh)                         ! Right half of target signature for hashing
                       Pt=LeftHash(S1)+RightHash(S2)              ! Storage position of target count
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime)   ! Add it up

                    CASE (Lower+4*Empty)
                       
                       ! Bottom edge is in state Lower and Top edge is in state Empty
                       ! Note that (E,L) and (L,E) has the same count after update by symmetry

                       ! Continue Lower arc end along bottom edge or continue Lower arc end along top edge

                       St=Sig-ISHFT(One,2*KL)+ISHFT(One,2*KU)  ! Target signature with Lower arc end at Top edge
                       S1=IAND(St,HashMask)                        
                       S2=ISHFT(St,-2*Wh)                         
                       Pt=LeftHash(S1)+RightHash(S2)               
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime) 
                       WalkCount(Ps)=WalkCount(Pt)



                    CASE (Upper+4*Empty)
                       
                       ! Bottom edge is in state Upper and Top edge is in state Empty
                       ! Note that (E,U) and (U,E) has the same count after update by symmetry

                       ! Continue Upper arc end along Bottom edge or continue Upper arc end along Top edge

                       St=Sig-ISHFT(One,2*KL+1)+ISHFT(One,2*KU+1)   ! Target signature with Upper arc end at Top edge
                       S1=IAND(St,HashMask)           
                       S2=ISHFT(St,-2*Wh)
                       Pt=LeftHash(S1)+RightHash(S2)
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime)
                       WalkCount(Ps)=WalkCount(Pt)

                    CASE (Lower+4*Lower)

                       ! Bottom edge is in state Lower and Top edge is in state Lower

                       ! Join Lower arc ends and relabel matching Upper arc end as Lower
                       St=Sig-5*ISHFT(One,2*KL)
                       CALL Relabel_Upper(St,KL)
                       S1=IAND(St,HashMask)
                       S2=ISHFT(St,-2*Wh)
                       Pt=LeftHash(S1)+RightHash(S2)
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime)



                    CASE (Upper+4*Lower)
                       
                       ! Bottom edge is in state Upper and Top edge is in state Lower

                       ! Join the two arc ends  
                       St=Sig-6*ISHFT(One,2*KL)
                       S1=IAND(St,HashMask)
                       S2=ISHFT(St,-2*Wh)
                       Pt=LeftHash(S1)+RightHash(S2)
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime)

                    CASE (Upper+4*Upper)
                       
                       ! Bottom edge is in state Upper and Top edge is in state Upper

                       ! Join Upper arc ends and relabel matching Lower arc end as Upper
                       St=Sig-10*ISHFT(One,2*KL)
                       CALL Relabel_Lower(St,KL)
                       S1=IAND(St,HashMask)
                       S2=ISHFT(St,-2*Wh)
                       Pt=LeftHash(S1)+RightHash(S2)
                       WalkCount(Pt)=MOD(WalkCount(Pt)+WalkCount(Ps),Prime)

                    END SELECT
                 END DO
              END IF
           END DO
        END DO

        !$OMP BARRIER

     END DO


     !$OMP BARRIER

     IF (Iam==0) THEN
        OPEN(7,FILE='hexsaw_cross_rhombus_omp.log',ACCESS='APPEND')
        CTime=OMP_GET_WTIME() 
        WRITE(7,'(A,I3,A,2F12.4)') "Calc -- Bot: ",Row,"  Wall Time: ",CTime-Ptime
        PTime=Ctime
        CLOSE(7)
        Ps=LeftHash(2)+RightHash(0)   ! Free end at bottom most edge
        IF (Row==Width) THEN
           WRITE(CW,'(I5)') Width
           SFN='hexsaw_cross_'//TRIM(ADJUSTL(CW))//'.data'
           K=LEN(SFN)
           OPEN(14,FILE=SFN(1:K),ACCESS='APPEND')
           WRITE(14,*) Prime,WalkCount(Ps)
           CLOSE(14)
        END IF
        CALL ClearSignatures
     END IF


  END DO
  !$OMP END PARALLEL



END SUBROUTINE ADDCOLUMNS




SUBROUTINE Relabel_Upper(S,Pos)

    ! Relabels an Upper arc end to Lower when joining two Upper arc ends

  USE saw_def
  IMPLICIT NONE
  INTEGER :: K,Sl,Pos,Num_Lowers
  INTEGER (KIND=8):: S

  Num_Lowers=0
  DO K=Pos,Width  ! Look for matching arc end above Pos
     Sl=IAND(ISHFT(S,-2*K),Three)
     IF (Sl==Lower) THEN
        Num_Lowers=Num_Lowers+1
     ELSE IF (Sl==Upper) THEN
        Num_Lowers=Num_Lowers-1
        IF (Num_Lowers<0) THEN  ! We found a match
           S=S-ISHFT(One,2*K)   ! now make it Lower
            RETURN
        END IF
     END IF
  END DO

  WRITE(7,*) 'Somethings wrong in Relabel_Upper ',Pos

END SUBROUTINE Relabel_Upper




SUBROUTINE Relabel_Lower(S,Pos)

  ! Relabels a Lower arc end to Upper when joining two Lower arc ends

  USE saw_def
  IMPLICIT NONE
  INTEGER :: K,Sl,Pos,Num_Uppers
  INTEGER (KIND=8):: S

  Num_Uppers=0
  DO K=Pos,0,-1  ! Look for matching Lower arc end below Pos
     Sl=IAND(ISHFT(S,-2*K),Three) 
     IF (Sl == Upper) THEN
        Num_Uppers=Num_Uppers+1
     ELSE IF (SL == Lower) THEN
        Num_Uppers=Num_Uppers-1
        IF (Num_Uppers < 0) THEN ! We found a match
           S=S+ISHFT(One,2*K)    ! now make it an Upper end
           RETURN
        END IF
     END IF
  END DO

  WRITE(7,*) 'Somethings wrong in Relabel_Lower',Pos
 
END SUBROUTINE Relabel_Lower


SUBROUTINE ConstructHashTable(Wl)

  USE saw_def
  IMPLICIT NONE
  INTEGER :: Wl,Wr,I,J,K,Signature(0:Width)
  INTEGER(KIND=8) :: Sig,Sl,Sr  ! Signatures stored as integers 

  Wr=Width+1-Wl

  ! Bases for packing all signatures (0,1)-->(Width,0) into a single array
  ! Just concatenating Motzkin paths (0,1)-->(Wl,I) with reversed Motzkin paths (0,0)-->(Wr,I)
  ! which gives all possible Motzkin paths (0,1)-->(Width,0).
  Base(0)=0
  DO I=0,MIN(Wl+1,Wr)
     Base(I+1)=Base(I)+Motzkin_1(Wl,I)*Motzkin_0(Wr,I) 
  END DO

  ALLOCATE(LeftHash(0:4**Wl))                ! Hash tables for Left and Right signatures
  ALLOCATE(RightHash(0:4**Wr))

  CALL ConstructSignatures(Wl)
 
  DO I=0,Wr  ! Height of Right signature
     K=0     ! Index number of signature at height I
     DO J=1,NumRightSignatures(I)  ! Loop over all signatures at height I
        K=K+1
        Sig=RightSignatures(J+BaseRightSignatures(I))  ! The coded signature 
        RightHash(Sig)=K       ! Value of Hash function is the index number
     END DO
  END DO

  DO I=0,Wl+1  ! Similar to above but for Left signatures
     K=0
     DO J=1,NumLeftSignatures(I)
        K=K+1
        Sig=LeftSignatures(J+BaseLeftSignatures(I))
        LeftHash(Sig)=Base(I)+(K-1)*Motzkin_0(Wr,I)  ! Value of Hash function for Left signature
     END DO
  END DO

  CALL ClearSignatures

END SUBROUTINE ConstructHashTable






SUBROUTINE ConstructThreadTable(Wl,Wm)

  USE saw_def
  IMPLICIT NONE
  INTEGER :: Wl,Wm,I,J,K,M,H,SigId
  INTEGER :: NumId(2**Wl),NumIndex(2**Wl),ThreadId(2**Wl),SignaturesPerThread(0:NumThreads-1)
  INTEGER(KIND=8) :: Sig,Sl,Sr  ! Signatures stored as integers

  TaskThread=0

  IF (Wl==Wm) THEN

     NumId=0
     DO H=0,Wl+1 

        DO J=1,NumLeftSignatures(H)
           Sig=LeftSignatures(J+BaseLeftSignatures(H))
           SigId=1
           DO K=0,Wm-1
              IF (IAND(Sig,ISHFT(Three,2*K))>0) SigId=SigId+ISHFT(1,K)
           END DO
           NumId(SigId)=NumId(SigId)+1

        END DO
     END DO
     CALL INDEXX(2**Wl,NumIndex,NumId)

     SignaturesPerThread=0
     ThreadId=-1
     I=2**Wl
     J=1
     K=0
     DO WHILE (I>=J)
        IF (NumId(NumIndex(I))>0) THEN
           ThreadId(NumIndex(I))=K
           SignaturesPerThread(K)=SignaturesPerThread(K)+NumId(NumIndex(I))
           DO WHILE (SignaturesPerThread(K)<SignaturesPerThread(0))
              IF (NumId(NumIndex(J))>0) THEN
                 ThreadId(NumIndex(J))=K
                 SignaturesPerThread(K)=SignaturesPerThread(K)+NumId(NumIndex(J))
              END IF
              J=J+1
              IF (I==J) EXIT
           END DO
           I=I-1
           K=MOD(K+1,NumThreads)
        ELSE
           I=I-1
        END IF
     END DO

     DO H=0,Wl+1  

        DO J=1,NumLeftSignatures(H)
           Sig=LeftSignatures(J+BaseLeftSignatures(H))   
           SigId=1
           DO K=0,Wm-1
              IF (IAND(Sig,ISHFT(Three,2*K))>0) SigId=SigId+ISHFT(1,K)
           END DO
           TaskThread(Sig)=ThreadId(SigId)
        END DO
     END DO

  ELSE

     NumId=0
     DO H=0,Wl+1 

        DO J=1,NumLeftSignatures(H)
           Sig=LeftSignatures(J+BaseLeftSignatures(H))
           SigId=1
           DO K=0,Wm-1
              IF (IAND(Sig,ISHFT(Three,2*K))>0) SigId=SigId+ISHFT(1,K)
           END DO
           DO K=Wm+2,Wl-1
              IF (IAND(Sig,ISHFT(Three,2*K))>0) SigId=SigId+ISHFT(1,K-2)
           END DO
           NumId(SigId)=NumId(SigId)+1

        END DO
     END DO

     CALL INDEXX(2**(Wl-2),NumIndex,NumId)

     SignaturesPerThread=0
     ThreadId=-1
     I=2**(Wl-2)
     J=1
     K=0
     DO WHILE (I>=J)
        IF (NumId(NumIndex(I))>0) THEN
           ThreadId(NumIndex(I))=K
           SignaturesPerThread(K)=SignaturesPerThread(K)+NumId(NumIndex(I))
           DO WHILE (SignaturesPerThread(K)<SignaturesPerThread(0))
              IF (NumId(NumIndex(J))>0) THEN
                 ThreadId(NumIndex(J))=K
                 SignaturesPerThread(K)=SignaturesPerThread(K)+NumId(NumIndex(J))
              END IF
              J=J+1
              IF (I==J) EXIT
           END DO
           I=I-1
           K=MOD(K+1,NumThreads)
        ELSE
           I=I-1
        END IF
     END DO

     DO H=0,Wl+1  
        DO J=1,NumLeftSignatures(H)
           Sig=LeftSignatures(J+BaseLeftSignatures(H))
           SigId=1
           DO K=0,Wm-1
              IF (IAND(Sig,ISHFT(Three,2*K))>0) SigId=SigId+ISHFT(1,K)
           END DO
           DO K=Wm+2,Wl-1
              IF (IAND(Sig,ISHFT(Three,2*K))>0) SigId=SigId+ISHFT(1,K-2)
           END DO
           TaskThread(Sig)=ThreadId(SigId)
        END DO
     END DO

  END IF



END SUBROUTINE ConstructThreadTable





SUBROUTINE ConstructSignatures(Wl)

  !     Construction of Left and Right signature as Motzkin paths
  !

  USE saw_def
  IMPLICIT NONE
  INTEGER, PARAMETER :: NumDir=3
  INTEGER :: K,Wl,Wr,Lat(0:Width), Path(0:Width),Steps, Pos 
  INTEGER :: NextDir(NumDir), CurDir
  INTEGER :: MaxSteps,I,K1,K2,M1,M2
  INTEGER, ALLOCATABLE :: Idx(:)
  INTEGER(KIND=8), ALLOCATABLE :: Signature(:)

  Wr=Width+1-Wl
  
  ! Bases for packing Left signatures into a single array
  BaseLeftSignatures(0)=0
  DO K=0,Wl+1
     BaseLeftSignatures(K+1)=BaseLeftSignatures(K)+Motzkin_1(Wl,K) ! There are Motzkin_1(Wl,K) paths (0,1)-->(Wl,K)
  END DO

  ! Bases for packing Right signatures into a single array
  BaseRightSignatures(0)=0
  DO K=0,Wr
     BaseRightSignatures(K+1)=BaseRightSignatures(K)+Motzkin_0(Wr,K) ! There are Motzkin_0(Wl,K) paths (0,0)-->(Wr,K)
  END DO

  ALLOCATE(LeftSignatures(BaseLeftSignatures(Wl+1)+1))   ! Array of Left signatures
  ALLOCATE(RightSignatures(BaseRightSignatures(Wr)+1))   ! Atrray of Right signaturea



  NextDir(1)=0
  NextDir(2)=1
  NextDir(3)=-1

  NumRightSignatures=0
  NumLeftSignatures=0

  ! Build the Left part of the signatures. Motzkin paths from (0,1) --> (MaxSteps,h)

  MaxSteps=Wl
  CALL BuildPaths(1,1)  ! Start with a Horizontal step (0,1) -->  (1,1)
  CALL BuildPaths(1,2)  ! Start with a Up step (0,1) -->  (1,2)
  CALL BuildPaths(1,3)  ! Start with a Down step (0,1) -->  (1,0)

  DO K2=0,MaxSteps+1       ! Take paths ending at height K2
     K1=NumLeftSignatures(K2)    ! there are K1 of these
     ALLOCATE(Idx(K1))     ! Create space for index and signature tables
     ALLOCATE(Signature(K1))
     M1=BaseLeftSignatures(K2)+1
     M2=BaseLeftSignatures(K2+1)
     Signature(1:K1)=LeftSignatures(M1:M2)  ! Here are the signatures
     CALL INDEXX8(K1,Idx,Signature)      ! now sort them in increasing order
     DO I=1,K1                      ! and put them back in the array of signatures
        LeftSignatures(M1+I-1)=Signature(Idx(I))
     END DO
     DEALLOCATE(Idx)              
     DEALLOCATE(Signature)
  END DO
  ! Left signatures now sorted in order for each height

  ! Build the Right part of the signatures. Reflected Motzkin paths from (0,0) --> (MaxSteps,h)

  MaxSteps=Wr
  CALL BuildPaths(0,1)  ! Start with a Horizontal step (0,0) -->  (1,0)
  CALL BuildPaths(0,2)  ! Start with a Up step (0,1) -->  (1,1)

  DO K2=0,MaxSteps  ! Same as above but for right paths
     K1=NumRightSignatures(K2)
     ALLOCATE(Idx(K1))
     ALLOCATE(Signature(K1))
     M1=BaseRightSignatures(K2)+1
     M2=BaseRightSignatures(K2+1)
     Signature(1:K1)=RightSignatures(M1:M2)
     CALL INDEXX8(K1,Idx,Signature)
     DO I=1,K1
        RightSignatures(M1+I-1)=Signature(Idx(I))
     END DO
     DEALLOCATE(Idx)
     DEALLOCATE(Signature)
  END DO
  ! Right signatures now sorted in order for each height
 

CONTAINS 


  SUBROUTINE BuildPaths(InitPos,InitStep)

    !    
    ! Build all Motzkin paths from (0,0) or (0,1) via basic back-tracking
    !     

    INTEGER :: M,InitPos,InitStep
    INTEGER(KIND=8), EXTERNAL :: LeftSig,RightSig
    LOGICAL :: Extended

    Lat=0
    Path=0
    Pos=InitPos
    Path(0)=Pos
    Lat(0)=NumDir
    Steps=1
    Pos=Pos+NextDir(InitStep)
    Path(Steps)=Pos
    DO WHILE (Steps>0) 
       DO WHILE (Steps<MaxSteps)
          CALL Extend_Path(Extended)
          IF (.NOT.Extended) THEN
             ! Contract Path
             Lat(Steps)=0
             Path(Steps)=0
             Steps=Steps-1
             Pos=Path(Steps)
             IF (Steps==0) RETURN
          END IF
       END DO
       IF (InitPos==1) THEN
          NumLeftSignatures(Pos)=NumLeftSignatures(Pos)+1
          M=BaseLeftSignatures(Pos)+NumLeftSignatures(Pos)
          LeftSignatures(M)=LeftSig(Path,Steps)
       ELSE
          NumRightSignatures(Pos)=NumRightSignatures(Pos)+1
          M=BaseRightSignatures(Pos)+NumRightSignatures(Pos)
          RightSignatures(M)=RightSig(Path,Steps)
       END IF

       ! Contract Path
       Lat(Steps)=0
       Path(Steps)=0
       Steps=Steps-1
       Pos=Path(Steps)
    END DO

  END SUBROUTINE BuildPaths


  SUBROUTINE Extend_Path(Extended)

    LOGICAL :: Extended

    CurDir=Lat(Steps)+1
    IF (CurDir>NumDir) THEN
       Extended=.FALSE.
    ELSE

       Pos=Pos+NextDir(CurDir)
       IF (Pos>=0) THEN
          Lat(Steps)=CurDir
          Steps=Steps+1
          Path(Steps)=Pos
          Extended=.TRUE.
       ELSE
          Extended=.FALSE.
       END IF

    END IF


  END SUBROUTINE Extend_Path




END SUBROUTINE ConstructSignatures



INTEGER(KIND=8) FUNCTION LeftSig(Path,Steps)

  USE saw_def
  IMPLICIT NONE
  INTEGER :: J, K, Steps,LocSignature(-1:1),Path(0:Width)

  ! Left signatures are Motzkin paths from height 1 

  LocSignature(1)=1
  LocSignature(0)=0
  LocSignature(-1)=2

  LeftSig=0
  K=1
  DO J=1,Steps
     LeftSig=LeftSig+K*LocSignature(Path(J)-Path(J-1))
     K=4*K
  END DO

END FUNCTION LeftSig


INTEGER(KIND=8) FUNCTION RightSig(Path,Steps)

  USE saw_def
  IMPLICIT NONE
  INTEGER :: J, K, Steps,LocSignature(-1:1),Path(0:Width)

  ! Right signatures are reflected Motzkin paths

  LocSignature(1)=2
  LocSignature(0)=0
  LocSignature(-1)=1

  RightSig=0
  K=1
  DO J=Steps,1,-1
     RightSig=RightSig+K*LocSignature(Path(J)-Path(J-1))
     K=4*K
  END DO

END FUNCTION RightSig



SUBROUTINE ClearSignatures

  USE saw_def
 
  DEALLOCATE(LeftSignatures)    
  DEALLOCATE(RightSignatures)    

END SUBROUTINE ClearSignatures



SUBROUTINE MotzkinNumbers()

  USE saw_def
  IMPLICIT NONE
  INTEGER :: I,J,M

  ! Number of Motzkin paths from (0,0) --> (I,J)

  Motzkin_0=0 
  Motzkin_0(0,0)=1
  DO I=0,Width+3
     DO J=0,I+1
        Motzkin_0(I+1,J)=Motzkin_0(I,J-1)+Motzkin_0(I,J)+Motzkin_0(I,J+1)
     END DO
  END DO

  ! Number of Motzkin paths from (0,1) --> (I,J)

  Motzkin_1=0
  Motzkin_1(0,1)=1
  DO I=0,Width+2
     DO J=0,I+2
        Motzkin_1(I+1,J)=Motzkin_1(I,J-1)+Motzkin_1(I,J)+Motzkin_1(I,J+1)
     END DO
  END DO



END SUBROUTINE MotzkinNumbers



SUBROUTINE PrintSignature(St,M,Mb)

  USE saw_def
  IMPLICIT NONE
  INTEGER :: I,J,M,Mb,Signature(Width)
  INTEGER(KIND=8) :: S,St,Sm
  
  Sm=Three
  if (Mb==2) Sm=One
  S=St
  DO I=1,M
     Signature(I)=IAND(S,Sm)
     S=S/Mb
  END DO
  WRITE(*,'(2X,30I1)') Signature(1:M)

END SUBROUTINE PrintSignature



SUBROUTINE INDEXX(MaxId,Idx,Signature)

  ! Sorting routine by indexing from Numerical Recipes 


  USE saw_def
  IMPLICIT NONE
  INTEGER ::  MaxId,Idx(MaxId),I,J,L,IR,INDXT
  INTEGER(KIND=4) ::  Signature(MaxId),Q

  DO J=1,MaxId
     Idx(J)=J
  END DO
  IF (MaxId==1) RETURN

  L=MaxId/2+1
  IR=MaxId
  DO WHILE (MaxId>0)
     IF (L>1) THEN
        L=L-1
        INDXT=Idx(L)
        Q=Signature(INDXT)
     ELSE
        INDXT=Idx(IR)
        Q=Signature(INDXT)
        Idx(IR)=Idx(1)
        IR=IR-1
        IF (IR==1) THEN
           Idx(1)=INDXT
           RETURN
        END IF
     END IF
     I=L
     J=L+L
     DO WHILE (J<=IR)
        IF (J<IR) THEN
           IF (Signature(Idx(J))<Signature(Idx(J+1))) J=J+1
        END IF
        IF (Q<Signature(Idx(J))) THEN
           Idx(I)=Idx(J)
           I=J
           J=J+J
        ELSE
           J=IR+1
        END IF
     END DO
     Idx(I)=INDXT
  END DO

END SUBROUTINE INDEXX




SUBROUTINE INDEXX8(MaxId,Idx,Signature)

  ! Sorting routine by indexing from Numerical Recipes 


  USE saw_def
  IMPLICIT NONE
  INTEGER ::  MaxId,Idx(MaxId),I,J,L,IR,INDXT
  INTEGER(KIND=8) ::  Signature(MaxId),Q

  DO J=1,MaxId
     Idx(J)=J
  END DO
  IF (MaxId==1) RETURN

  L=MaxId/2+1
  IR=MaxId
  DO WHILE (MaxId>0)
     IF (L>1) THEN
        L=L-1
        INDXT=Idx(L)
        Q=Signature(INDXT)
     ELSE
        INDXT=Idx(IR)
        Q=Signature(INDXT)
        Idx(IR)=Idx(1)
        IR=IR-1
        IF (IR==1) THEN
           Idx(1)=INDXT
           RETURN
        END IF
     END IF
     I=L
     J=L+L
     DO WHILE (J<=IR)
        IF (J<IR) THEN
           IF (Signature(Idx(J))<Signature(Idx(J+1))) J=J+1
        END IF
        IF (Q<Signature(Idx(J))) THEN
           Idx(I)=Idx(J)
           I=J
           J=J+J
        ELSE
           J=IR+1
        END IF
     END DO
     Idx(I)=INDXT
  END DO

END SUBROUTINE INDEXX8

