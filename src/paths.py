from __future__ import print_function

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils'))
from files_and_paths import Urls

#Host = Urls.metadataWebService
Host = "http://192.168.1.46:9008/"
BaseWwwTmpDir = os.path.join(os.path.dirname(__file__), '../www-tmp')
BaseWwwDir = os.path.join(os.path.dirname(__file__), '../www')
