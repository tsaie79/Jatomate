from fireworks import explicit_serialize, FiretaskBase, FWAction

from atomate.utils.utils import env_chk, get_logger, logger
from atomate.vasp.database import VaspCalcDb

from pyzfs.common.wfc.vasploader import VaspWavefunctionLoader
from pyzfs.zfs.main import ZFSCalculation

from pymatgen.io.vasp.inputs import Structure

from monty.serialization import loadfn
from monty.json import jsanitize

import subprocess
import os
import json


@explicit_serialize
class RunPyzfs(FiretaskBase):
    """
    run pyzfs
    zfs_cmd:
        srun -n 2048 -c 2 python /home/tug03990/site-packages/pyzfs/examples/VASP/run.py > out (cori)
        mpiexec -n 20 pyzfs --wfcfmt vasp > out (owls, efrc)
    """
    required_params = ["pyzfs_cmd"]

    def run_task(self, fw_spec):

        wd = os.getcwd()

        try:
            raw_struct = Structure.from_file(wd + "/POSCAR")
            formula = raw_struct.composition.formula
            structure = raw_struct.as_dict()

        except:
            formula = None
            structure = None

        cmd = env_chk(self["pyzfs_cmd"], fw_spec)
        logger.info("Running command: {}".format(cmd))
        return_code = subprocess.call([cmd], shell=True)
        logger.info("Command {} finished running with returncode: {}".format(cmd, return_code))

        return FWAction(
            update_spec={
                "structure": structure,
                "formula": formula,
            }
        )


@explicit_serialize
class PyzfsToDb(FiretaskBase):

    optional_params = ["db_file", "additional_fields", "collection_name"]

    def run_task(self, fw_spec):

        pyzfs_out = loadfn("pyzfs_out.json")
        pyzfs_out = jsanitize(pyzfs_out)

        additional_fields = self.get("additional_fields", {})
        d = additional_fields.copy()
        d["formula"] = fw_spec["formula"]
        d["structure"] = fw_spec["structure"]
        d["pyzfs_out"] = pyzfs_out
        d["dir_name"] = os.getcwd()
        # store the results
        db_file = env_chk(self.get("db_file"), fw_spec)
        if not db_file:
            with open("pyzfs_todb.json", "w") as f:
                f.write(json.dumps(d, default=DATETIME_HANDLER, indent=4))
        else:
            db = VaspCalcDb.from_db_file(db_file, admin=True)
            print(self.get("collection_name", db.collection.name))
            db.collection = db.db[self.get("collection_name", db.collection.name)]
            t_id = db.insert(d)
            logger.info("Pyzfs calculation complete.")
        return FWAction()
