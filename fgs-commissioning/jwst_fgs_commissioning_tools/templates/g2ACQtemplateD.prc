WAIT ((ILRV(IFGS_CMD_ACC_CNT) != accCount) || (ILRV(IFGS_CMD_REJ_CNT) != rejCount))

IF (ILRV(IFGS_CMD_REJ_CNT) != rejCount)
   SYSALERT INFO, "[TEST] @IFGS_CONFIG was rejected for Guider 2 (ACQ2)"
   PRINTF "\n@IFGS_CONFIG was rejected for Guider 2 (ACQ2)\n"
   PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
   SUSPEND
ELSE
    WAIT (ELRV(IFGS_CONFIG_STAT) != 'STARTED')

    IF (ELRV(IFGS_CONFIG_STAT) != 'SUCCESS')
       SYSALERT INFO, "[TEST] @IFGS_CONFIG did not complete successfully for Guider 2 (ACQ2)"
       PRINTF "\n@IFGS_CONFIG did not complete successfully for Guider 2 (ACQ2)\n"
       PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
       SUSPEND
    ENDIF
ENDIF


;------------------------------------------------------------------------------
; Image Generation - Guider 2
;------------------------------------------------------------------------------

; Wait for FGSES readiness

PRINTF "Press GO to expose, or ABORT ALL to abort the procedure\n"
SUSPEND

; Begin exposure

SYSALERT INFO, "[TEST] Beginning exposure on Guider 2..."

accCount = ILRV(IFGS_CMD_ACC_CNT)
rejCount = ILRV(IFGS_CMD_REJ_CNT)

@IFGS_EXPOSE GUIDER2

WAIT ((ILRV(IFGS_CMD_ACC_CNT) != accCount) || (ILRV(IFGS_CMD_REJ_CNT) != rejCount))

IF (ILRV(IFGS_CMD_REJ_CNT) != rejCount)
   SYSALERT INFO, "[TEST] @IFGS_EXPOSE was rejected for Guider 2"
   PRINTF "\n@IFGS_EXPOSE was rejected for Guider 2\n"
   PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
   SUSPEND
ELSE
    WAIT (ELRV(IFGS_EXPOSE_STAT) != 'STARTED')

    IF (ELRV(IFGS_EXPOSE_STAT) != 'SUCCESS')
       SYSALERT INFO, "[TEST] @IFGS_EXPOSE did not complete successfully for Guider 2"
       PRINTF "\n@IFGS_EXPOSE did not complete successfully for Guider 2\n"
       PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
       SUSPEND
    ENDIF
ENDIF

SYSALERT INFO, "[TEST] ACQing on Guider 2..."
PRINTF "\nACQing on Guider 2...\n"
PRINTF "To return the FSW to STANDBY on Guider 2, press Go\n"
PRINTF "To abort the procedure, press ABORT ALL\n"
SUSPEND


;------------------------------------------------------------------------------
; Return to STANDBY - Guider 2
;------------------------------------------------------------------------------

; Move to STANDBY

SYSALERT INFO, "[TEST] Moving to STANDBY on Guider 2..."

accCount = ILRV(IFGS_CMD_ACC_CNT)
rejCount = ILRV(IFGS_CMD_REJ_CNT)

@IFGS_TRANSITION GUIDER2, STANDBY

WAIT ((ILRV(IFGS_CMD_ACC_CNT) != accCount) || (ILRV(IFGS_CMD_REJ_CNT) != rejCount))

IF (ILRV(IFGS_CMD_REJ_CNT) != rejCount)
   SYSALERT INFO, "[TEST] @IFGS_TRANSITION was rejected for Guider 2 (to STANDBY)"
   PRINTF "\n@IFGS_TRANSITION was rejected for Guider 2 (to STANDBY)\n"
   PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
   SUSPEND
ELSE
    WAIT (ELRV(IFGS_TRANS_STAT) != 'STARTED')

    IF (ELRV(IFGS_TRANS_STAT) != 'SUCCESS')
       SYSALERT INFO, "[TEST] @IFGS_TRANSITION did not complete successfully for Guider 2 (to STANDBY)"
       PRINTF "\n@IFGS_TRANSITION did not complete successfully for Guider 2 (to STANDBY)\n"
       PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
       SUSPEND
    ENDIF
ENDIF

; Abort exposures and free FPAPs

@IDAQ_ABORT 1
WAIT (1)

@IDAQ_ABORT 2
WAIT (1)

@IDAQ_ABORT 3
WAIT (1)

	; Send the data into the pipeline
	;st gisim_FVTS_EXPOSE_COMPLETE

;------------------------------------------------------------------------------
; Procedure complete
;------------------------------------------------------------------------------

SYSALERT INFO, "[TEST] ACQ G2 testing complete"
PRINTF "\nACQ G2 testing complete\n"

ENDP
