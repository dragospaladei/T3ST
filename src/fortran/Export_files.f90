MODULE output_files

  USE constants, ONLY: rp

  IMPLICIT NONE
  PRIVATE

! ---- named file indices (sequential, no renumbering needed) ----
INTEGER, PARAMETER, PUBLIC :: f_parameters = 1

! ---- initial state ----
INTEGER, PARAMETER, PUBLIC :: f_Xinit      = f_parameters + 1
INTEGER, PARAMETER, PUBLIC :: f_Yinit      = f_Xinit      + 1
INTEGER, PARAMETER, PUBLIC :: f_Zinit      = f_Yinit      + 1
INTEGER, PARAMETER, PUBLIC :: f_Vpinit     = f_Zinit      + 1
INTEGER, PARAMETER, PUBLIC :: f_muinit     = f_Vpinit     + 1
INTEGER, PARAMETER, PUBLIC :: f_betaFinit  = f_muinit     + 1
INTEGER, PARAMETER, PUBLIC :: f_betaDinit  = f_betaFinit  + 1
INTEGER, PARAMETER, PUBLIC :: f_FMinit     = f_betaDinit  + 1

! ---- trajectories ----
INTEGER, PARAMETER, PUBLIC :: f_Xtraj      = f_FMinit     + 1
INTEGER, PARAMETER, PUBLIC :: f_Ytraj      = f_Xtraj      + 1
INTEGER, PARAMETER, PUBLIC :: f_Ztraj      = f_Ytraj      + 1
INTEGER, PARAMETER, PUBLIC :: f_Vptraj     = f_Ztraj      + 1
INTEGER, PARAMETER, PUBLIC :: f_mutraj     = f_Vptraj     + 1
INTEGER, PARAMETER, PUBLIC :: f_betaFtraj  = f_mutraj     + 1
INTEGER, PARAMETER, PUBLIC :: f_betaDtraj  = f_betaFtraj  + 1
INTEGER, PARAMETER, PUBLIC :: f_FMtraj     = f_betaDtraj  + 1
INTEGER, PARAMETER, PUBLIC :: f_Htraj      = f_FMtraj     + 1
INTEGER, PARAMETER, PUBLIC :: f_q1traj     = f_Htraj      + 1
INTEGER, PARAMETER, PUBLIC :: f_q2traj     = f_q1traj     + 1
INTEGER, PARAMETER, PUBLIC :: f_q3traj     = f_q2traj     + 1
INTEGER, PARAMETER, PUBLIC :: f_Pctraj     = f_q3traj     + 1
INTEGER, PARAMETER, PUBLIC :: f_ck1traj    = f_Pctraj     + 1
INTEGER, PARAMETER, PUBLIC :: f_ck2traj    = f_ck1traj    + 1
INTEGER, PARAMETER, PUBLIC :: f_ck3traj    = f_ck2traj    + 1

! ---- transport ----
INTEGER, PARAMETER, PUBLIC :: f_Obs_GF     = f_ck3traj    + 1
INTEGER, PARAMETER, PUBLIC :: f_Obs_FF     = f_Obs_GF     + 1
INTEGER, PARAMETER, PUBLIC :: f_Obs_FD     = f_Obs_FF     + 1
INTEGER, PARAMETER, PUBLIC :: f_Obs_DF     = f_Obs_FD     + 1
INTEGER, PARAMETER, PUBLIC :: f_Obs_DD     = f_Obs_DF     + 1

! ---- final state ----
INTEGER, PARAMETER, PUBLIC :: f_Xout       = f_Obs_DD     + 1
INTEGER, PARAMETER, PUBLIC :: f_Yout       = f_Xout       + 1
INTEGER, PARAMETER, PUBLIC :: f_Zout       = f_Yout       + 1
INTEGER, PARAMETER, PUBLIC :: f_Vpout      = f_Zout       + 1
INTEGER, PARAMETER, PUBLIC :: f_muout      = f_Vpout      + 1
INTEGER, PARAMETER, PUBLIC :: f_betaFout   = f_muout      + 1
INTEGER, PARAMETER, PUBLIC :: f_betaDout   = f_betaFout   + 1
INTEGER, PARAMETER, PUBLIC :: f_FMout      = f_betaDout   + 1
INTEGER, PARAMETER, PUBLIC :: f_Hout       = f_FMout      + 1
INTEGER, PARAMETER, PUBLIC :: f_q1out      = f_Hout       + 1
INTEGER, PARAMETER, PUBLIC :: f_q2out      = f_q1out      + 1

! ---- correlations ----
INTEGER, PARAMETER, PUBLIC :: f_Vcorff     = f_q2out      + 1
INTEGER, PARAMETER, PUBLIC :: f_VcorTT     = f_Vcorff     + 1
INTEGER, PARAMETER, PUBLIC :: f_VcorTN     = f_VcorTT     + 1
INTEGER, PARAMETER, PUBLIC :: f_VcorNT     = f_VcorTN     + 1

INTEGER, PARAMETER, PUBLIC :: f_Lagr_corr  = f_VcorNT     + 1
INTEGER, PARAMETER, PUBLIC :: f_mask       = f_Lagr_corr  + 1

! ---- total number of files ----
INTEGER, PARAMETER, PUBLIC :: nfiles = f_mask
  ! ---- public file units ----
  INTEGER, PUBLIC :: funit(nfiles)

  ! ---- internal bookkeeping ----
  CHARACTER(len=128) :: fname(nfiles)
  CHARACTER(len=8)   :: fext(nfiles)
  LOGICAL            :: enabled(nfiles)
  LOGICAL            :: is_open(nfiles)

  PUBLIC :: init_output_files
  PUBLIC :: open_output_files
  PUBLIC :: close_output_files
  PUBLIC :: enable_output
  PUBLIC :: disable_output

  PUBLIC :: write_initial_state
  PUBLIC :: write_trajectories
  PUBLIC :: write_transport
  PUBLIC :: write_final_state
  PUBLIC :: write_correlations
  PUBLIC :: write_lagrangian
  PUBLIC :: write_mask
CONTAINS

  SUBROUTINE init_output_files()
    INTEGER :: i

    DO i = 1, nfiles
       funit(i)   = 100 + i
       fname(i)   = ''
       fext(i)    = '.bin'
       enabled(i) = .TRUE.
       is_open(i) = .FALSE.
    END DO

    fname(f_parameters) = 'parameters'
    fext(f_parameters)  = '.dat'

    fname(f_Xinit)      = 'Xinit'
    fname(f_Yinit)      = 'Yinit'
    fname(f_Zinit)      = 'Zinit'
    fname(f_Vpinit)     = 'Vpinit'
    fname(f_muinit)     = 'muinit'
    fname(f_betaFinit)  = 'betaFinit'
    fname(f_betaDinit)  = 'betaDinit'
    fname(f_FMinit)     = 'FMinit'

    fname(f_Xtraj)      = 'Xtraj'
    fname(f_Ytraj)      = 'Ytraj'
    fname(f_Ztraj)      = 'Ztraj'
    fname(f_Vptraj)     = 'Vptraj'
    fname(f_mutraj)     = 'mutraj'
    fname(f_betaFtraj)  = 'betaFtraj'
    fname(f_betaDtraj)  = 'betaDtraj'
    fname(f_FMtraj)     = 'FMtraj'
    fname(f_Htraj)      = 'Htraj'
    fname(f_q1traj)     = 'q1traj'
    fname(f_q2traj)     = 'q2traj'
    fname(f_q3traj)     = 'q3traj'
    fname(f_Pctraj)     = 'Pctraj'
    fname(f_ck1traj)    = 'ck1traj'
    fname(f_ck2traj)    = 'ck2traj'
    fname(f_ck3traj)    = 'ck3traj'

    fname(f_Obs_GF)     = 'Obs_GF'
    fname(f_Obs_FF)     = 'Obs_FF'
    fname(f_Obs_FD)     = 'Obs_FD'
    fname(f_Obs_DF)     = 'Obs_DF'
    fname(f_Obs_DD)     = 'Obs_DD'

    fname(f_Xout)       = 'Xout'
    fname(f_Yout)       = 'Yout'
    fname(f_Zout)       = 'Zout'
    fname(f_Vpout)      = 'Vpout'
    fname(f_muout)      = 'muout'
    fname(f_betaFout)   = 'betaFout'
    fname(f_betaDout)   = 'betaDout'
    fname(f_FMout)      = 'FMout'
    fname(f_Hout)       = 'Hout'
    fname(f_q1out)      = 'q1out'
    fname(f_q2out)      = 'q2out'

    fname(f_Vcorff)     = 'Vcorff'
    fname(f_VcorTT)     = 'VcorTT'
    fname(f_VcorTN)     = 'VcorTN'
    fname(f_VcorNT)     = 'VcorNT'
    fname(f_Lagr_corr)  = 'Lagr_corr'
    fname(f_mask)       = 'mask'

  END SUBROUTINE init_output_files


  SUBROUTINE open_output_files(folder)
    CHARACTER(len=*), INTENT(IN) :: folder

    INTEGER :: i, ios
    CHARACTER(len=:), ALLOCATABLE :: output_folder
    CHARACTER(len=:), ALLOCATABLE :: path

    output_folder = TRIM(folder)
    IF (LEN_TRIM(output_folder) == 0) THEN
       WRITE(*,*) 'ERROR: empty output folder'
       STOP
    END IF

    IF (output_folder(LEN_TRIM(output_folder):LEN_TRIM(output_folder)) /= '/') THEN
       output_folder = output_folder//'/'
    END IF

    DO i = 1, nfiles

       IF (.NOT. enabled(i)) CYCLE
       IF (LEN_TRIM(fname(i)) == 0) CYCLE

       path = output_folder//TRIM(fname(i))//TRIM(fext(i))

       IF (TRIM(fext(i)) == '.dat') THEN
          OPEN(UNIT=funit(i), FILE=TRIM(path), &
               STATUS='replace', ACTION='write', IOSTAT=ios)
       ELSE
          OPEN(UNIT=funit(i), FILE=TRIM(path), &
               FORM='unformatted', ACCESS='stream', &
               STATUS='replace', ACTION='write', IOSTAT=ios)
       END IF

       IF (ios /= 0) THEN
          WRITE(*,*) 'ERROR: could not open output file: ', TRIM(path)
          WRITE(*,*) 'IOSTAT = ', ios
          STOP
       END IF

       is_open(i) = .TRUE.

    END DO

  END SUBROUTINE open_output_files


  SUBROUTINE close_output_files()
    INTEGER :: i, ios

    DO i = 1, nfiles

       IF (is_open(i)) THEN
          CLOSE(UNIT=funit(i), IOSTAT=ios)

          IF (ios /= 0) THEN
             WRITE(*,*) 'WARNING: could not close file unit ', funit(i)
             WRITE(*,*) 'IOSTAT = ', ios
          END IF

          is_open(i) = .FALSE.
       END IF

    END DO

  END SUBROUTINE close_output_files


  SUBROUTINE enable_output(file_id)
    INTEGER, INTENT(IN) :: file_id

    IF (.NOT. valid_file_id(file_id)) THEN
       WRITE(*,*) 'ERROR: invalid file id in enable_output: ', file_id
       STOP
    END IF

    enabled(file_id) = .TRUE.

  END SUBROUTINE enable_output


  SUBROUTINE disable_output(file_id)
    INTEGER, INTENT(IN) :: file_id

    IF (.NOT. valid_file_id(file_id)) THEN
       WRITE(*,*) 'ERROR: invalid file id in disable_output: ', file_id
       STOP
    END IF

    enabled(file_id) = .FALSE.

  END SUBROUTINE disable_output


  LOGICAL FUNCTION valid_file_id(file_id)
    INTEGER, INTENT(IN) :: file_id

    valid_file_id = file_id >= 1 .AND. file_id <= nfiles

  END FUNCTION valid_file_id

SUBROUTINE write_initial_state(X, Y, Z, Vp, mu, betaF, betaD, FM)
   REAL(rp), INTENT(IN) :: X(:), Y(:), Z(:)
   REAL(rp), INTENT(IN) :: Vp(:), mu(:), betaF(:), betaD(:), FM(:)

   WRITE(funit(f_Xinit))  X
   WRITE(funit(f_Yinit))  Y
   WRITE(funit(f_Zinit))  Z
   WRITE(funit(f_Vpinit)) Vp
   WRITE(funit(f_muinit)) mu
   WRITE(funit(f_betaFinit)) betaF
   WRITE(funit(f_betaDinit)) betaD
   WRITE(funit(f_FMinit))    FM
END SUBROUTINE write_initial_state


SUBROUTINE write_trajectories(Xtraj, Ytraj, Ztraj, Vptraj, mutraj, betaFtraj, betaDtraj, FMtraj, Htraj, &
                              q1traj, q2traj, q3traj, Pctraj, ck1traj, ck2traj, ck3traj)
   REAL(rp), INTENT(IN) :: Xtraj(:, :), Ytraj(:, :), Ztraj(:, :)
   REAL(rp), INTENT(IN) :: Vptraj(:, :), mutraj(:, :)
   REAL(rp), INTENT(IN) :: betaFtraj(:, :), betaDtraj(:, :), FMtraj(:, :)
   REAL(rp), INTENT(IN) :: Htraj(:, :), Pctraj(:, :)
   REAL(rp), INTENT(IN) :: q1traj(:, :), q2traj(:, :), q3traj(:, :)
   REAL(rp), INTENT(IN) :: ck1traj(:, :), ck2traj(:, :), ck3traj(:, :)

   WRITE(funit(f_Xtraj))   Xtraj
   WRITE(funit(f_Ytraj))   Ytraj
   WRITE(funit(f_Ztraj))   Ztraj
   WRITE(funit(f_Vptraj))  Vptraj
   WRITE(funit(f_mutraj))  mutraj
   WRITE(funit(f_betaFtraj)) betaFtraj
   WRITE(funit(f_betaDtraj)) betaDtraj
   WRITE(funit(f_FMtraj))    FMtraj
   WRITE(funit(f_Htraj))   Htraj

   WRITE(funit(f_q1traj))  q1traj
   WRITE(funit(f_q2traj))  q2traj
   WRITE(funit(f_q3traj))  q3traj

   WRITE(funit(f_Pctraj))  Pctraj
   WRITE(funit(f_ck1traj)) ck1traj
   WRITE(funit(f_ck2traj)) ck2traj
   WRITE(funit(f_ck3traj)) ck3traj
END SUBROUTINE write_trajectories


SUBROUTINE write_transport(Obs_GF, Obs_FF, Obs_FD, Obs_DF, Obs_DD)
   REAL(rp), INTENT(IN) :: Obs_GF(:, :), Obs_FF(:, :), Obs_FD(:, :), Obs_DF(:, :), Obs_DD(:, :)

   WRITE(funit(f_Obs_GF)) Obs_GF
   WRITE(funit(f_Obs_FF)) Obs_FF
   WRITE(funit(f_Obs_FD)) Obs_FD
   WRITE(funit(f_Obs_DF)) Obs_DF
   WRITE(funit(f_Obs_DD)) Obs_DD
END SUBROUTINE write_transport


SUBROUTINE write_final_state(X, Y, Z, Vp, mu, betaF, betaD, FM, Ham, q1al, q2al)
   REAL(rp), INTENT(IN) :: X(:), Y(:), Z(:)
   REAL(rp), INTENT(IN) :: Vp(:), mu(:), betaF(:), betaD(:), FM(:), Ham(:)
   REAL(rp), INTENT(IN) :: q1al(:), q2al(:)

   WRITE(funit(f_Xout))  X
   WRITE(funit(f_Yout))  Y
   WRITE(funit(f_Zout))  Z
   WRITE(funit(f_Vpout)) Vp
   WRITE(funit(f_muout)) mu
   WRITE(funit(f_betaFout)) betaF
   WRITE(funit(f_betaDout)) betaD
   WRITE(funit(f_FMout))    FM
   WRITE(funit(f_Hout))  Ham
   WRITE(funit(f_q1out)) q1al
   WRITE(funit(f_q2out)) q2al
END SUBROUTINE write_final_state


SUBROUTINE write_correlations(Vcorff, VcorTT, VcorTN, VcorNT)
   REAL(rp), INTENT(IN) :: Vcorff(:, :), VcorTT(:, :)
   REAL(rp), INTENT(IN) :: VcorTN(:, :), VcorNT(:, :)

   WRITE(funit(f_Vcorff)) Vcorff
   WRITE(funit(f_VcorTT)) VcorTT
   WRITE(funit(f_VcorTN)) VcorTN
   WRITE(funit(f_VcorNT)) VcorNT
END SUBROUTINE write_correlations


SUBROUTINE write_lagrangian(Lagr_corr)
   REAL(rp), INTENT(IN) :: Lagr_corr(:)

   WRITE(funit(f_Lagr_corr)) Lagr_corr
END SUBROUTINE write_lagrangian


SUBROUTINE write_mask(mask_value)
   REAL(rp), INTENT(IN) :: mask_value

   WRITE(funit(f_mask)) mask_value
END SUBROUTINE write_mask


END MODULE output_files
