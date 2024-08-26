

      SUBROUTINE CHEM_FIELD( CGRID, JDATE, JTIME, FND_SPC )

C-----------------------------------------------------------------------
C Function: Read a species field from an input file and replace current
C           CGRID values. Used for OH in a reduced chemistry model.
C Revision History:
C     16 Aug 2024 S.L.Napelenok: initial implementation
C-----------------------------------------------------------------------              

      USE GRID_CONF, ONLY:  NCOLS, NROWS, NLAYS
      USE CGRID_SPCS, ONLY: N_CGRID_SPC, CGRID_NAME
!     USE RUNTIME_VARS, ONLY: LOGDEV, OUTDEV, CMYPE, APPL_NAME
      USE RUNTIME_VARS, ONLY: OUTDEV, CMYPE, APPL_NAME
      USE UTILIO_DEFN
!     USE CENTRALIZED_IO_MODULE, ONLY: interpolate_var
      USE CENTRALIZED_IO_MODULE, ONLY: STRTCOLSTD, ENDCOLSTD, STRTROWSTD, ENDROWSTD
#ifndef mpas
#ifdef parallel
      USE SE_MODULES          ! stenex (using SE_GLOBAL_SUM_MODULE)
#else
      USE NOOP_MODULES        ! stenex (using NOOP_GLOBAL_SUM_MODULE)
#endif
#endif

      IMPLICIT NONE

C...Arguments:
      REAL, POINTER :: CGRID( :,:,:,: )  ! Species concentrations
      INTEGER JDATE           ! Current date (YYYYDDD)
      INTEGER JTIME           ! Current time (HHMMSS)
      CHARACTER( 16 ) FND_SPC ! Species to replace

C...Local Variables:
      CHARACTER( 16 ), SAVE :: PNAME = 'CHEM_FIELD'
      LOGICAL, SAVE :: FIRSTIME = .TRUE. ! flag for first pass thru
      INTEGER     ASTAT
      CHARACTER( 120 ) :: XMSG = ' '
      REAL( 4 ), ALLOCATABLE, SAVE :: CHEM_FIELD_ARRAY ( :,:,: )
      INTEGER LAY, ROW, COL
      INTEGER, SAVE :: SPC_IND = -1
      CHARACTER(16), SAVE :: CHEM_FIELD_FILE = 'CHEM_FIELD_FILE'

C-----------------------------------------------------------------------

C...Initialization
      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.  

         ALLOCATE ( CHEM_FIELD_ARRAY ( NCOLS,NROWS,NLAYS ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating CHEM_FIELD_ARRAY'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         SPC_IND = INDEX1( FND_SPC, N_CGRID_SPC, CGRID_NAME )

         IF ( SPC_IND .LT. 1 ) THEN
            XMSG = 'Failure finding requested species' // TRIM(FND_SPC)
     &             // 'in CGRID array. (CHEM_FIELD)'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

C...Open input file
         IF ( .NOT. OPEN3( CHEM_FIELD_FILE, FSREAD3, PNAME ) ) THEN
            XMSG = 'Could not open CHEM_FIELD_FILE'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         ENDIF

         IF ( .NOT. DESC3( CHEM_FIELD_FILE ) ) THEN
            XMSG = 'Could not get ' //CHEM_FIELD_FILE// 
     &             ' file description'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      END IF

C-----------------------------------------------------------------------

C...Read data
      CHEM_FIELD_ARRAY = 0.0

!     CALL INTERPOLATE_VAR (FND_SPC, JDATE, JTIME, CHEM_FIELD_ARRAY, CHEM_FIELD_FILE )

      IF ( .NOT. INTERPX( CHEM_FIELD_FILE, FND_SPC, PNAME,
     &                    STRTCOLSTD,  ENDCOLSTD, STRTROWSTD, ENDROWSTD, 1, NLAYS,
     &                    JDATE, JTIME, CHEM_FIELD_ARRAY ) ) THEN
         XMSG = ' Error interpolating variable from ' // CHEM_FIELD_FILE
         Call M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF

      DO LAY = 1, NLAYS
         DO ROW = 1, NROWS
            DO COL = 1, NCOLS
               CGRID( COL,ROW,LAY,SPC_IND ) = 
     &                   CHEM_FIELD_ARRAY ( COL,ROW,LAY )
            END DO
         END DO

      END DO

      RETURN
      END
              
