
WAIT ((ILRV(IFGS_CMD_ACC_CNT) != accCount) || (ILRV(IFGS_CMD_REJ_CNT) != rejCount))

IF (ILRV(IFGS_CMD_REJ_CNT) != rejCount)
   SYSALERT INFO, "[TEST] @IFGS_CONFIG was rejected for Guider 1 (ACQ1)"
   PRINTF "\n@IFGS_CONFIG was rejected for Guider 1 (ACQ1)\n"
   PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
   SUSPEND
ELSE
    WAIT (ELRV(IFGS_CONFIG_STAT) != 'STARTED')

    IF (ELRV(IFGS_CONFIG_STAT) != 'SUCCESS')
       SYSALERT INFO, "[TEST] @IFGS_CONFIG did not complete successfully for Guider 1 (ACQ1)"
       PRINTF "\n@IFGS_CONFIG did not complete successfully for Guider 1 (ACQ1)\n"
       PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
       SUSPEND
    ENDIF
ENDIF

; Configure FPE for ACQ2

accCount = ILRV(IFGS_CMD_ACC_CNT)
rejCount = ILRV(IFGS_CMD_REJ_CNT)

