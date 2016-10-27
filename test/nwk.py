#!/usr/bin/env python

import sys
import re

nwk = open(sys.argv[1], 'r').read().strip()
print nwk
quoted = re.sub(r'([^();,]+)', '"\\1"', nwk.replace(";", ""))
print quoted
print eval(quoted)

