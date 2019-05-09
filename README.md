# HST_LIGO_DD


prob_trigger.py [--options]
returns the likelihood that we should execute a trigger
options are:
   --verbose (shows all output, by default)
   --quiet (shows only the likelihood)
   --help
   --plotit (generates plots)
   --area = [float] (in deg2)
   --distance = [float] (in Mpc)
   --bprob = [float] (BNS probablity)
   --velocity = [float] (spectral expansion velocity, if known)
   --hours = [float] (time since burst, in hours)
Options can be combined. There are some default values for testing
