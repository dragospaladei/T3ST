
    PROGRAM DNS
    IMPLICIT NONE

!========================================================================================================
!============ VARIABLE DECLARATION ======================================================================
!========================================================================================================
    CHARACTER(LEN=10000)                         :: CWD                             ! current working directory (for exporting data)
    INTEGER                                      :: lenCWD                          ! length of CWD
    CHARACTER(LEN=:), ALLOCATABLE                :: filein,fileout                  ! address where data is exported
    CHARACTER(LEN=35)                            :: time, shot                      ! auxiliary, for folder definition
    INTEGER                                      :: i, j, nlines, iostat,ngrid1,ngrid2! index of the simulation
    CHARACTER(LEN=:), ALLOCATABLE                :: folder                          ! name of the folder for exporting data
    INTEGER                                      :: run ! index of the simulation

    REAL            , ALLOCATABLE, DIMENSION(:,:):: Vx     ! positions, parallel velocity, magnetic moment, probability - for RK4
    REAL            , ALLOCATABLE, DIMENSION(:)  :: z     ! positions, parallel velocity, magnetic moment, probability - for RK4
    CHARACTER(LEN=80)                            :: zk                                 ! current working directory (for exporting data)
!========================================================================================================
!============ ADDRESS DECLARATION FOR IMPORT DATA ====================================================
!========================================================================================================
    CALL GetCWD(CWD)
    filein = TRIM(ADJUSTL(CWD))
    fileout = filein
    lenCWD = LEN(filein)
    shot = '54719'
    shot = '52702'
    shot = '54178'
    
    !============ run = time of data ====================================================
!do run = 1,15  ! for 54719
do run = 1,15   ! for 52702
    write(time, '(I0)') run
!============ file to import and export ====================================================
    filein = filein(1:lenCWD-13)//'\'//trim(shot)//'\geqdsk_'//trim(shot)//'_'//trim(time)//'_201x201'      
    if(run.lt.10) then
        fileout = fileout(1:lenCWD-13)//'\'//trim(shot)//'\geqdsk_'//trim(shot)//'_00'//trim(time)//'.dat'        
    elseif((10.le.run).and.(run.lt.100)) then
        fileout = fileout(1:lenCWD-13)//'\'//trim(shot)//'\geqdsk_'//trim(shot)//'_0'//trim(time)//'.dat'        
    elseif(100.le.run) then
        fileout = fileout(1:lenCWD-13)//'\'//trim(shot)//'\geqdsk_'//trim(shot)//'_'//trim(time)//'.dat'        
    endif
    iostat = 0
    open(file=trim(filein), status='old', action='read', iostat=iostat, unit=100)
    open(file=trim(fileout), action='write', iostat=iostat, unit=200)
!============ file dimension ====================================================
    write(*,*) filein
    write(*,*) fileout
!============ read first line to extract grid dimension ====================================================
    read(100, '(A)', iostat=iostat) zk
    read(zk(55:56), '(I)', iostat=iostat) ngrid1
    read(zk(59:61), '(I)', iostat=iostat) ngrid2
!============ see how many lines are there; the first is lost from the previous import ====================================================
    nlines = 0
    do while (iostat == 0)
            read(100, '(A)', iostat=iostat) zk
            nlines = nlines + 1
    enddo
!============ start from beggining to import ====================================================
    rewind(unit=100)

    write(*,*) 'here', nlines
!============ remove the first line ====================================================
    read(100, '(A)') zk

    do i=1,nlines-2      
      read(100, '(A)') zk
      do j = 0,4 
         write(200,*) zk(16*j+1:16*(j+1))
      enddo
    enddo
close(100)
close(200)
write(*,*) run
enddo
    END PROGRAM
