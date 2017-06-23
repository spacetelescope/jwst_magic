
WAIT ((ILRV(IFGS_CMD_ACC_CNT) != accCount) || (ILRV(IFGS_CMD_REJ_CNT) != rejCount))

IF (ILRV(IFGS_CMD_REJ_CNT) != rejCount)
   SYSALERT INFO, "[TEST] @IFGS_GUIDESTAR was rejected for Guider 1"
   PRINTF "\n@IFGS_GUIDESTAR was rejected for Guider 1"
   PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
   SUSPEND
ELSE
    WAIT (ELRV(IFGS_GUIDESTAR_STAT) != 'STARTED')

    IF (ELRV(IFGS_GUIDESTAR_STAT) != 'SUCCESS')
       SYSALERT INFO, "[TEST] @IFGS_GUIDESTAR did not complete successfully for Guider 1"
       PRINTF "\n@IFGS_GUIDESTAR did not complete successfully for Guider 1\n"
       PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
       SUSPEND
    ENDIF
ENDIF

	; Set up the pipeline and the SSR for the data
	;st gisim_spw_clrcnts
	;st gisim_FVTS_EXPOSE_SETUP('FVTS_FGSES_ACQ_G1')


;------------------------------------------------------------------------------
; DAQ configuration
;------------------------------------------------------------------------------

	; **** THIS CONFIG WORKS ****    
	@IDAQ_CONFIGMSG OBSID="FVTS_ACQ1_G1",      \
	                EXPID=6,                   \
	                NINTS=3,                   \
	                AVERAGE=1,                 \
	                NGROUPS=2,                 \
	                NFRAMES=1,                 \
	                NROWS=128,                 \
	                NCOLS=128,                 \
	                FRAME1=DROP,               \
	                DETECTOR=17,               \
	                SLOT=1,                    \
	                DATAMODE=1,                \
	                RETAIN=STORE,              \
	                NIMAGES=1,                 \
	                NOTIFY=2,                  \
	                NROWSPKT=1,                \
                        BITS2SHIFT = 0
                     ;  NINTS_AVG=1
	
	; Wait for data acq to ackowldege that it is ready
	wait  (ELRV(IDAQ_EXP_STATUS1) == 'SETUP_COMPLETE')
	
        wait 10
        spaceWireAddr1 = ILRV(IDAQ_GUIDER1_ADDR)
        groupNum1 = ILRV(IDAQ_GUIDER1_GPI)

	; Tell DAQ to save the file
	@IDAQ_SAVE  EXTENSION='FIM',    \
	            FILE='ACQ1RRRR',        \
	            FILESTORE='fgs:/',  \
	            SCEP=10
	wait 2

   	@IDAQ_CONFIGMSG OBSID="FVTS_ACQ2_G1",   \
                   EXPID=6,                \
                   NINTS=5,                \
                   NGROUPS=2,              \
                   NFRAMES=1,              \
                   NROWS=32,               \
                   NCOLS=32,               \
                   FRAME1=DROP,            \
                   DETECTOR=17,            \
                   SLOT=2,                 \
                   DATAMODE=1,             \
                   RETAIN=STORE,           \
                   NIMAGES=1,              \
                   NOTIFY=2,               \
                   NROWSPKT=4,             \
                   AVERAGE=1,              \
                   BITS2SHIFT=0
		 ; NINTS_AVG=1

    ; Wait for data acq to ackowldege that it is ready
    wait  (ELRV(IDAQ_EXP_STATUS2) == 'SETUP_COMPLETE')

        wait 10
        spaceWireAddr2 = ILRV(IDAQ_GUIDER1_ADDR)
        groupNum2 = ILRV(IDAQ_GUIDER1_GPI)

    ; Tell DAQ to save the file
    @IDAQ_SAVE  EXTENSION='FIM', FILE='ACQ2CDS', FILESTORE='fgs:/', SCEP=11
    wait 2


; Configure FPE for ACQ1

accCount = ILRV(IFGS_CMD_ACC_CNT)
rejCount = ILRV(IFGS_CMD_REJ_CNT)

