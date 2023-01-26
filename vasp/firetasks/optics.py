import subprocess

from pydash.objects import has, get

from fireworks import FiretaskBase, FWAction, explicit_serialize
from fireworks.utilities.fw_serializers import DATETIME_HANDLER

from pymatgen.io.vasp.inputs import *
from pymatgen.io.vasp.sets import MPStaticSet, MVLGWSet, MPHSEBSSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

from atomate.vasp.database import VaspCalcDb
from atomate.utils.utils import env_chk
from atomate.vasp.config import *
from atomate.vasp.drones import VaspDrone
from atomate.common.firetasks.glue_tasks import get_calc_loc



from monty.shutil import compress_dir, decompress_dir

from glob import glob

import shutil, gzip, os, re, traceback, time


@explicit_serialize
class WriteMVLGWFromPrev(FiretaskBase):
    """
    Writes input files for a static run. Assumes that output files from a
    previous (e.g., optimization) run can be accessed in current dir or
    prev_calc_dir. Also allows lepsilon (dielectric constant) calcs.
    Optional params:
        potcar_spec (bool): Instead of writing the POTCAR, write a
            "POTCAR.spec". This is intended to allow testing of workflows
            without requiring pseudo-potentials to be installed on the system.
        (documentation for all other optional params can be found in
        MPStaticSet)
    """

    optional_params = [
        "prev_calc_dir",
        "prev_incar",
        "nbands",
        "reciprocal_density",
        "mode",
        "nbands_factor",
        "ncores",
        "other_params"
    ]

    def run_task(self, fw_spec):
        other_params = self.get("other_params", {})
        vis = MVLGWSet.from_prev_calc(
            prev_calc_dir=self.get("prev_calc_dir", "."),
            prev_incar=self.get("prev_incar", None),
            nbands=self.get("nbands", None),
            reciprocal_density=self.get("reciprocal_density", 100),
            mode=self.get("mode", "DIAG"),
            copy_wavecar=False,
            nbands_factor=self.get("nbands_factor", 5),
            ncores=self.get("ncores", 16),
            **other_params
        )

        vis.write_input(".")