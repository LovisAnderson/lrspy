import logging
DETAILED_DEBUG_LEVEL_NUM = 9
logging.addLevelName(DETAILED_DEBUG_LEVEL_NUM, 'DetailedDebug')

def detailed_debug(self, message, *args, **kws):
    if self.isEnabledFor(DETAILED_DEBUG_LEVEL_NUM):
        # Yes, logger takes its '*args' as 'args'.
        self._log(DETAILED_DEBUG_LEVEL_NUM, message, args, **kws)

logging.Logger.detailed_debug = detailed_debug