
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

; Configure reference stars

accCount = ILRV(IFGS_CMD_ACC_CNT)
rejCount = ILRV(IFGS_CMD_REJ_CNT)

