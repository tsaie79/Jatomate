from fireworks import explicit_serialize, FiretaskBase, FWAction

from atomate.utils.utils import env_chk, get_logger, logger
from atomate.vasp.database import VaspCalcDb

from pymatgen.io.vasp.inputs import Structure

from monty.serialization import loadfn
from monty.json import jsanitize

import subprocess
import os
import json
from pydash.objects import has, get


@explicit_serialize
class RunPyzfs(FiretaskBase):
    """
    For run zfs, check github "pyzfs/examples/VASP/"
    zfs_cmd:
        srun -n 2048 -c 2 python ~/site-packages/pyzfs/examples/VASP/run.py > out (cori) # of node = 64
        mpiexec -n 100 pyzfs --wfcfmt vasp > out (owls, efrc) #!! -n 100 is needed!!
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

    optional_params = ["db_file", "additional_fields", "collection_name", "task_fields_to_push"]

    def run_task(self, fw_spec):

        pyzfs_out = loadfn("pyzfs_out.json")
        pyzfs_out = jsanitize(pyzfs_out)

        additional_fields = self.get("additional_fields", {})
        d = additional_fields.copy()
        d["formula"] = fw_spec["formula"]
        d["structure"] = fw_spec["structure"]
        d["pyzfs_out"] = pyzfs_out
        d["dir_name"] = os.getcwd()
        # Automatically add prev fws information
        for prev_info_key in ["prev_fw_taskid", "prev_fw_db", "prev_fw_collection"]:
            if prev_info_key in fw_spec:
                d.update({prev_info_key: fw_spec[prev_info_key]})

        # store the results
        db_file = env_chk(self.get("db_file"), fw_spec)
        if not db_file:
            with open("pyzfs_todb.json", "w") as f:
                f.write(json.dumps(d, default=DATETIME_HANDLER, indent=4))
        else:
            db = VaspCalcDb.from_db_file(db_file, admin=True)
            print(self.get("collection_name", db.collection.name))
            db.collection = db.db[self.get("collection_name", db.collection.name)]
            d.update({"db": db.db_name, "collection": db.collection.name})
            t_id = db.insert(d)
            logger.info("Pyzfs calculation complete.")

        return FWAction()
