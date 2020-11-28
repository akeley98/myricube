#!/usr/bin/env python3
# Convert stdout log of makefile to format suitable for rtags's rc.
# Usage: make --dry-run | cckiss/rtags.py | rc -c -
import sys, os
for line in sys.stdin:
    if line.startswith("cckiss/cckiss"):
        # See cckiss_rtags_main function of cckiss.cc .
        os.system("CCKISS_RTAGS_PRINT=1 %s" % line)
    else:
        sys.stdout.write(line)
