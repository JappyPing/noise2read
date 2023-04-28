# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-04-28 10:05:05
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-04-28 10:21:22
#!/usr/bin/env python
import sys

if sys.version_info[0] < 3 or (sys.version_info[0] >= 3 and sys.version_info[1] < 5):
    raise RuntimeError('Must be using Python 3.5 or later')

import SourceLibrary.IsoMiRmap_Library as START
if __name__ == "__main__":
    START.start_generation()
