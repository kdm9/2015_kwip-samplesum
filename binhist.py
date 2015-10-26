#!/usr/bin/env python
from __future__ import print_function, division
from collections import Counter
import khmer
import sys

LEVELS = [0, 1, 2, 10, 128, 254, 255]
print('File', '\t'.join(str(i) for i in LEVELS), sep='\t')
for ctgz in sys.argv[1:]:
    ct = khmer.load_countgraph(ctgz)
    table = ct.get_raw_tables()[0]
    ctr = Counter()
    for count in table:
        ctr[ord(count)] +=1
    hist = {i:0 for i in LEVELS}
    for count, freq in ctr.items():
        for level in LEVELS:
            if count <= level:
                hist[level] += freq
                break
    print(ctgz, '\t'.join(str(v) for k, v in sorted(hist.items())), sep='\t')
