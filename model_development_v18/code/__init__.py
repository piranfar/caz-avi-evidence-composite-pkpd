"""hujam-caz-avi — uncertainty-integrated joint target attainment for
continuous-infusion ceftazidime-avibactam, with a decision-analytic layer
(value of information, misselection/regret, triage) built to be portable to
other beta-lactam/beta-lactamase-inhibitor pairs.

Installed under the import name ``hujam`` (see ../pyproject.toml); this
directory is also usable directly as a script folder without installation,
which is how it was developed and tested.

Start here:
    ../SOFTWARE.md            what this is, install, quickstart, how to adapt
    interface.py              the two-function engine contract + conformance check
    model2_engine.py          the CAZ-AVI-specific engine (verified against v16)
    model2_hujam.py           the portable decision layer (VOI, regret, EVPPI)
    joint_popk_nlme.py        Model 1: the clearance-correlation estimate
    test_model1.py            119 checks covering Model 1 and the primary model
    test_model2.py            decision-layer identities and sampler checks
"""

__version__ = "1.0.0"

# The modules in this package import their siblings with bare names
# (`import model2_engine as E`) rather than relative imports, because the
# package was developed and is still primarily used as a directory of scripts
# run directly (`python model2_hujam.py`), not as an installed library. That
# only works if this directory is on sys.path. Running a script directly
# guarantees that automatically; `import hujam.something` after installation
# does not, so it is done explicitly here, once, on package import.
import os as _os
import sys as _sys

_here = _os.path.dirname(_os.path.abspath(__file__))
if _here not in _sys.path:
    _sys.path.insert(0, _here)
