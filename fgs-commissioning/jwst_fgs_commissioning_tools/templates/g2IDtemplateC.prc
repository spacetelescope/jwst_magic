WAIT ((ILRV(IFGS_CMD_ACC_CNT) != accCount) || (ILRV(IFGS_CMD_REJ_CNT) != rejCount))

IF (ILRV(IFGS_CMD_REJ_CNT) != rejCount)
   SYSALERT INFO, "[TEST] @IFGS_REFCOUNT was rejected for Guider 2"
   PRINTF "\n@IFGS_REFCOUNT was rejected for Guider 2"
   PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
   SUSPEND
ELSE
    WAIT (ELRV(IFGS_REFCNT_STAT) != 'STARTED')

    IF (ELRV(IFGS_REFCNT_STAT) != 'SUCCESS')
       SYSALERT INFO, "[TEST] @IFGS_REFCOUNT did not complete successfully for Guider 2"
       PRINTF "\n@IFGS_REFCOUNT did not complete successfully for Guider 2\n"
       PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
       SUSPEND
    ENDIF
ENDIF

accCount = ILRV(IFGS_CMD_ACC_CNT)
rejCount = ILRV(IFGS_CMD_REJ_CNT)

