#!/usr/bin/env python3

import sys
import pandas as pd
print(sys.argv)

df = pd.read_csv(sys.argv[1], sep="\t")
df.astype(int).to_csv(sys.argv[2], sep="\t")