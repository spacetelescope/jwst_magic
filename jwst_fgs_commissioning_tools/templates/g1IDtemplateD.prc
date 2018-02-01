
WAIT ((ILRV(IFGS_CMD_ACC_CNT) != accCount) || (ILRV(IFGS_CMD_REJ_CNT) != rejCount))

IF (ILRV(IFGS_CMD_REJ_CNT) != rejCount)
   SYSALERT INFO, "[TEST] @IFGS_REFSTAR was rejected for Guider 1"
   PRINTF "\n@IFGS_REFSTAR was rejected for Guider 1"
   PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
   SUSPEND
ELSE
	wait 2
ENDIF

accCount = ILRV(IFGS_CMD_ACC_CNT)
rejCount = ILRV(IFGS_CMD_REJ_CNT)

