
WAIT ((ILRV(IFGS_CMD_ACC_CNT) != accCount) || (ILRV(IFGS_CMD_REJ_CNT) != rejCount))

IF (ILRV(IFGS_CMD_REJ_CNT) != rejCount)
   SYSALERT INFO, "[TEST] @IFGS_REFSTAR was rejected for Guider 1"
   PRINTF "\n@IFGS_REFSTAR was rejected for Guider 1"
   PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
   SUSPEND
ELSE
    WAIT (ELRV(IFGS_REFSTAR_STAT) != 'STARTED')

    IF (ELRV(IFGS_REFSTAR_STAT) != 'SUCCESS')
       SYSALERT INFO, "[TEST] @IFGS_REFSTAR did not complete successfully for Guider 1"
       PRINTF "\n@IFGS_REFSTAR did not complete successfully for Guider 1\n"
       PRINTF "Press GO to continue, or ABORT ALL to abort the procedure\n"
       SUSPEND
    ENDIF
ENDIF

	; Set up the pipeline and the SSR for the data
	;st gisim_spw_clrcnts
	;st gisim_FVTS_EXPOSE_SETUP('FVTS_IDDATA_FGSES')


