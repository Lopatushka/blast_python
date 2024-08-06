import os
import subprocess
import io
import numpy as np
import pandas as pd
import warnings
from collections import Counter
from Bio import AlignIO
from Bio import SeqIO

s = "./db/16S_rRNAs_fecal_long/16S_rRNAs_fecal_long"
print(s.split("/")[-1])
