import numpy as np
import pandas as pd
import os

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
#  pip install rpy2
import argparse


# do we actually want to run the python file from here tho?
# only works if we are only using the python code and not the R code? --> make it compatible?
