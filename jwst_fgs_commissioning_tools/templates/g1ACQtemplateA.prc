

;------------------------------------------------------------------------------
; Variable Declaration
;------------------------------------------------------------------------------

STRING msgText   ; Will store the text to display in the RTL

INT accCount   ; Holds the value of the accepted counter
INT rejCount   ; Holds the value of the rejected counter

; Config parameters
INT spaceWireAddr1
INT spaceWireAddr2
INT groupNum1
INT groupNum2

SYSALERT INFO, "[TEST] Beginning FGSES ACQ test..."

	;ST GISIM_FVTS_SIDESEL('p','p')  
	;ST FVTS_PLT_FPESETUP('p')
	@ICFM_FORMAT "fgs:/"
	wait 2

;------------------------------------------------------------------------------
; FPE setup
;------------------------------------------------------------------------------

; Put the FSW into FPE Operation state
if (elrv("IFGS_ENGG1_ICESTAT") != "ICE_IDLE")
  @IFGS_ICEOP GUIDER1
endif

WAIT 10

if (elrv("IFGS_ENGG1_FPESTATE") != "FPE_SIDECAR_OPERATION")
  @IFGS_ASICOP GUIDER1
endif

WAIT 10


; Reset to STANDBY

accCount = ILRV(IFGS_CMD_ACC_CNT)
rejCount = ILRV(IFGS_CMD_REJ_CNT)

@IFGS_TRANSITION GUIDER1, STANDBY

WAIT ((ILRV(IFGS_CMD_ACC_CNT) != accCount) || (ILRV(IFGS_CMD_REJ_CNT) != rejCount))

IF (ILRV(IFGS_CMD_REJ_CNT) != rejCount)
   SYSALERT INFO, "[TEST] @IFGS_TRANSITION was rejected for Guider 1 (STANDBY)"
   PRINTF "@IFGS_TRANSITION was rejected for Guider 1 (STANDBY)\n"
   PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
   SUSPEND
ELSE
    WAIT (ELRV(IFGS_TRANS_STAT) != 'STARTED')

    IF (ELRV(IFGS_TRANS_STAT) != 'SUCCESS')
       SYSALERT INFO, "[TEST] @IFGS_TRANSITION did not complete successfully for Guider 1 (STANDBY)"
       PRINTF "@IFGS_TRANSITION did not complete successfully for Guider 1 (STANDBY)\n"
       PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
       SUSPEND
    ENDIF
ENDIF

; Transition to ACQ

accCount = ILRV(IFGS_CMD_ACC_CNT)
rejCount = ILRV(IFGS_CMD_REJ_CNT)

@IFGS_TRANSITION GUIDER1, ACQ

WAIT ((ILRV(IFGS_CMD_ACC_CNT) != accCount) || (ILRV(IFGS_CMD_REJ_CNT) != rejCount))

IF (ILRV(IFGS_CMD_REJ_CNT) != rejCount)
   SYSALERT INFO, "[TEST] @IFGS_TRANSITION was rejected for Guider 1 (to ACQ)"
   PRINTF "\n@IFGS_TRANSITION was rejected for Guider 1 (to ACQ)\n"
   PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
   SUSPEND
ELSE
    WAIT (ELRV(IFGS_TRANS_STAT) != 'STARTED')

    IF (ELRV(IFGS_TRANS_STAT) != 'SUCCESS')
       SYSALERT INFO, "[TEST] @IFGS_TRANSITION did not complete successfully for Guider 1 (to ACQ)"
       PRINTF "\n@IFGS_TRANSITION did not complete successfully for Guider 1 (to ACQ)\n"
       PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
       SUSPEND
    ENDIF
ENDIF

;------------------------------------------------------------------------------
; FPE configuration - Guider 1
;------------------------------------------------------------------------------

SYSALERT INFO, "[TEST] Configuring Guider 1..."

; Configure guide star

accCount = ILRV(IFGS_CMD_ACC_CNT)
rejCount = ILRV(IFGS_CMD_REJ_CNT)

