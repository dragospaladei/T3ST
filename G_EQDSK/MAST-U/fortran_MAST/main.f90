
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
    shot = '47090'
    shot = '47000'
    shot = '47020'
    shot = '47980'
    shot = '47089'
    shot = '49200'
!============ run = time of data ====================================================
do run = 5,49 
    write(time, '(I0)') run*20
!============ file to import and export ====================================================
    filein = filein(1:lenCWD-13)//'\'//trim(shot)//'\g0'//trim(shot)//'.00'//trim(time)                                          ! must have a 'source' and a 'data' folder (in the same folder as source)
    fileout = fileout(1:lenCWD-13)//'\'//trim(shot)//'\g0'//trim(shot)//'_'//trim(time)//'.dat'                                          ! must have a 'source' and a 'data' folder (in the same folder as source)
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
write(*,*) run*20
enddo
    END PROGRAM
