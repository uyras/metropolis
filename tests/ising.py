import subprocess
import sys
from io import StringIO
import numpy as np


result = subprocess.run(
    sys.argv[2:-1],
    cwd=sys.argv[1],
    capture_output=True,
    text=True,
    timeout=5
)
assert result.returncode == 0

dat_orig = np.loadtxt(StringIO(result.stdout))
dat = np.loadtxt(sys.argv[1]+sys.argv[-1])
assert np.isclose(dat[:,[0,1,2]], dat_orig[:,[0,1,2]]).all()

assert result.stderr == ""