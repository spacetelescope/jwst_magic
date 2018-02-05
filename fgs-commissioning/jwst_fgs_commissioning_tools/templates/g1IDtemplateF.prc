; Set up the pipeline and the SSR for the data
;st gisim_spw_clrcnts
;st gisim_FVTS_EXPOSE_SETUP('FVTS_IDDATA_FGSES')

;------------------------------------------------------------------------------
; DAQ configuration
;------------------------------------------------------------------------------

	; **** THIS CONFIG WORKS ****    
	print 'FVTS_IDDATA: Configure data acquisition for Guider 1...'
	@IDAQ_CONFIGMSG    OBSID="FVTS_IDDATA_G1", \
	            EXPID=6,                       \
	            NINTS=72,                      \
	            AVERAGE=1,                     \
	            NGROUPS=2,                     \
	            NFRAMES=1,                     \
	            NROWS=64,                      \
	            NCOLS=2048,                    \
	            FRAME1=DROP,                   \
	            DETECTOR=17,                   \
	            SLOT=1,                        \
	            DATAMODE=1,                    \
	            RETAIN=STORE,                  \
	            NIMAGES=2,                     \
	            NOTIFY=2,                      \
	            NROWSPKT=1,                    \
             	    BITS2SHIFT = 0
                    NINTS_AVG=1

	                      
	; Wait for data acq to ackowldege that it is ready
	wait  (ELRV(IDAQ_EXP_STATUS1) == 'SETUP_COMPLETE')

        wait 10
        spaceWireAddr = ILRV(IDAQ_GUIDER1_ADDR)
        groupNum = ILRV(IDAQ_GUIDER1_GPI)

	; Tell DAQ to save the file
	@IDAQ_SAVE  EXTENSION='FIM', \
	    FILE='ID1STRIP',         \
	    FILESTORE='fgs:/',       \
	    SCEP=10

	wait 2

	; Configure the SWTS to flow the science data
	 @gisim_SPW2_CMD CMD='swts clrcnts'
	wait 1


; Configure FPE for ID

accCount = ILRV(IFGS_CMD_ACC_CNT)
rejCount = ILRV(IFGS_CMD_REJ_CNT)

@IFGS_CONFIG DETECTOR=GUIDER1,        \
             SWADDRESS=spaceWireAddr, \
             SLOT=1,                  \
             NINTS=1,                 \
             NGROUPS=groupNum,        \
             NFRAMES=1,               \
             NSAMPLES=1,              \
             GROUPGAP=1,              \
             NROWS=64,                \
             NCOLS=2048,              \
             ROWCORNER=-70.5113,      \
             COLCORNER=70.8061
          ;  NINTS_AVG=1

WAIT ((ILRV(IFGS_CMD_ACC_CNT) != accCount) || (ILRV(IFGS_CMD_REJ_CNT) != rejCount))

IF (ILRV(IFGS_CMD_REJ_CNT) != rejCount)
   SYSALERT INFO, "[TEST] @IFGS_CONFIG was rejected for Guider 2 (ID)"
   PRINTF "\n@IFGS_CONFIG was rejected for Guider 2 (ID)\n"
   PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
   SUSPEND
ELSE
    WAIT (ELRV(IFGS_CONFIG_STAT) != 'STARTED')

    IF (ELRV(IFGS_CONFIG_STAT) != 'SUCCESS')
       SYSALERT INFO, "[TEST] @IFGS_CONFIG did not complete successfully for Guider 2 (ID)"
       PRINTF "\n@IFGS_CONFIG did not complete successfully for Guider 2 (ID)\n"
       PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
       SUSPEND
    ENDIF
ENDIF


;------------------------------------------------------------------------------
; Image Generation - Guider 1
;------------------------------------------------------------------------------

; Wait for FGSES readiness

PRINTF "Press GO to expose, or ABORT ALL to abort the procedure\n"
SUSPEND

; Begin exposure

SYSALERT INFO, "[TEST] Beginning exposure on Guider 1..."

accCount = ILRV(IFGS_CMD_ACC_CNT)
rejCount = ILRV(IFGS_CMD_REJ_CNT)

@IFGS_EXPOSE GUIDER1

WAIT ((ILRV(IFGS_CMD_ACC_CNT) != accCount) || (ILRV(IFGS_CMD_REJ_CNT) != rejCount))

IF (ILRV(IFGS_CMD_REJ_CNT) != rejCount)
   SYSALERT INFO, "[TEST] @IFGS_EXPOSE was rejected for Guider 1"
   PRINTF "\n@IFGS_EXPOSE was rejected for Guider 1\n"
   PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
   SUSPEND
ELSE
    WAIT (ELRV(IFGS_EXPOSE_STAT) != 'STARTED')

    IF (ELRV(IFGS_EXPOSE_STAT) != 'SUCCESS')
       SYSALERT INFO, "[TEST] @IFGS_EXPOSE did not complete successfully for Guider 1"
       PRINTF "\n@IFGS_EXPOSE did not complete successfully for Guider 1\n"
       PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
       SUSPEND
    ENDIF
ENDIF

SYSALERT INFO, "[TEST] IDing on Guider 1..."
PRINTF "\nIDing on Guider 1...\n"
PRINTF "To return the FSW to STANDBY on Guider 1, press Go\n"
PRINTF "To abort the procedure, press ABORT ALL\n"
SUSPEND


;------------------------------------------------------------------------------
; Return to STANDBY - Guider 1
;------------------------------------------------------------------------------

; Move to STANDBY

SYSALERT INFO, "[TEST] Moving to STANDBY on Guider 1..."

accCount = ILRV(IFGS_CMD_ACC_CNT)
rejCount = ILRV(IFGS_CMD_REJ_CNT)

@IFGS_TRANSITION GUIDER1, STANDBY

WAIT ((ILRV(IFGS_CMD_ACC_CNT) != accCount) || (ILRV(IFGS_CMD_REJ_CNT) != rejCount))

IF (ILRV(IFGS_CMD_REJ_CNT) != rejCount)
   SYSALERT INFO, "[TEST] @IFGS_TRANSITION was rejected for Guider 1 (to STANDBY)"
   PRINTF "\n@IFGS_TRANSITION was rejected for Guider 1 (to STANDBY)\n"
   PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
   SUSPEND
ELSE
    WAIT (ELRV(IFGS_TRANS_STAT) != 'STARTED')

    IF (ELRV(IFGS_TRANS_STAT) != 'SUCCESS')
       SYSALERT INFO, "[TEST] @IFGS_TRANSITION did not complete successfully for Guider 1 (to STANDBY)"
       PRINTF "\n@IFGS_TRANSITION did not complete successfully for Guider 1 (to STANDBY)\n"
       PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
       SUSPEND
    ENDIF
ENDIF

; Abort exposures and free FPAPs

@IDAQ_RELEASE 1
wait 2

@IDAQ_ABORT 1
WAIT (1)

	; Send the data into the pipeline
	;st gisim_FVTS_EXPOSE_COMPLETE
	wait 5

;------------------------------------------------------------------------------
; Procedure complete
;------------------------------------------------------------------------------

SYSALERT INFO, "[TEST] ID G1 testing complete"
PRINTF "\nID G1 testing complete\n"

ENDP
