# Guide for doing some code stuffs

### Start logging
I recommend using logging for all outputs that are not strictly for debugging, but even still
it can be nice to be able to look at output logs. There is a pre prepped system for you to use
a preconfigured logger.

``` python
import logging
from src.General.setenv import setupLogging
setupLogging()
```

This will allow you to log outputs with 

``` python
logging.info(<string>)
logging.debug(<string>)
logging.warning(<string>)
```

and it will be preserved in a log file while also popping up on stdout
also for Winston, I made some icecream versions of it lol. You can look
in the code for them :)

### icecream
This is a good way to do debugging, but shouldn't be kept in for the final product. Super easy to use

``` python
from icecream import ic
ic(<var>)
```

totally bonkers. use the `ic(<var>)` instead of print.

### download any new required modules
just run makeenv.py, easy as that