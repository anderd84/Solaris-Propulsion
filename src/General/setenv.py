import logging
import datetime
import os
import sys
from icecream import ic

def setupLogging():
    if not logging.getLogger().hasHandlers():
        time = datetime.datetime.now()
        os.makedirs("./logs", exist_ok=True)

        logger = logging.getLogger()
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter("[%(levelname)s] : %(message)s")

        fileH = logging.FileHandler(filename=f"./logs/{time.strftime('[%d-%m-%Y] %H-%M-%S')}.log")
        fileH.setFormatter(formatter)

        console = logging.StreamHandler(sys.stdout)
        console.setFormatter(formatter)

        logger.addHandler(console)
        logger.addHandler(fileH)

        logging.debug("Logger initialized")
    else:
        logging.warning("Logger already initialized")

def icDebug(debugOutput):
    logging.debug(ic.format(debugOutput))

def icInfo(infoOutput):
    logging.debug(ic.format(infoOutput))