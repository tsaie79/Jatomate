"""
Firetasks for FWs.

"""

import shutil
import json
import os
from pydash.objects import has, get
import numpy as np

from monty.json import MontyEncoder, jsanitize

from spglib import standardize_cell

from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Incar, Outcar, Kpoints

from pytopomat.irvsp_caller import IRVSPCaller, IRVSPOutput, IRVSPOutputAll
from pytopomat.vasp2trace_caller import (
    Vasp2TraceCaller,
    Vasp2Trace2Caller,
    Vasp2TraceOutput,
)
from pytopomat.z2pack_caller import Z2PackCaller

from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.utilities.fw_serializers import DATETIME_HANDLER

from atomate.utils.utils import env_chk, get_logger
from atomate.utils.database import CalcDb
from atomate.vasp.database import VaspCalcDb


logger = get_logger(__name__)


@explicit_serialize
class RunIRVSP(FiretaskBase):
    """
    Execute IRVSP in current directory.

    """
    optional_params = ["set_spn", "symprec"]
    def run_task(self, fw_spec):

        wd = os.getcwd()
        set_spn = self["set_spn"]
        symprec = self["symprec"]
        irvsp_caller = IRVSPCaller(wd, set_spn=set_spn, symprec=symprec)

        try:
            raw_struct = Structure.from_file(wd + "/POSCAR")
            formula = raw_struct.composition.formula
            structure = raw_struct.as_dict()

            outcar = Outcar(wd + "/OUTCAR")
            efermi = outcar.efermi

        except:
            formula = None
            structure = None
            efermi = None

        kpoints = Kpoints.from_file(wd + "/KPOINTS")
        data = IRVSPOutput(wd + "/outir.txt", kpoints)

        return FWAction(
            update_spec={
                "irvsp_out": data.as_dict(),
                "structure": structure,
                "formula": formula,
                "efermi": efermi,
                "post_relax_sg_name": irvsp_caller.sg_name,
                "post_relax_sg_number": irvsp_caller.sg_number
            }
        )

@explicit_serialize
class RunIRVSPAll(FiretaskBase):
    """
    Execute IRVSP in current directory.

    """
    required_params = ["set_spn", "symprec"]
    def run_task(self, fw_spec):

        wd = os.getcwd()
        set_spn = self["set_spn"]
        symprec = self["symprec"]
        irvsp_caller = IRVSPCaller(wd, set_spn=set_spn, symprec=symprec)

        try:
            raw_struct = Structure.from_file(wd + "/POSCAR")
            formula = raw_struct.composition.formula
            structure = raw_struct.as_dict()

            outcar = Outcar(wd + "/OUTCAR")
            efermi = outcar.efermi

        except:
            formula = None
            structure = None
            efermi = None

        kpoints = Kpoints.from_file(wd + "/KPOINTS")
        general = IRVSPOutputAll(wd + "/outir.txt")
        high_sym_data = IRVSPOutput(wd + "/outir.txt", kpoints)
        data = general.as_dict().copy()
        data["parity_eigenvals"] = {"high_sym": high_sym_data.parity_eigenvals, "general": general.parity_eigenvals}

        return FWAction(
            update_spec={
                "irvsp_out": data,
                "structure": structure,
                "formula": formula,
                "efermi": efermi,
                "post_relax_sg_name": irvsp_caller.sg_name,
                "post_relax_sg_number": irvsp_caller.sg_number
            }
        )

@explicit_serialize
class RunIRVSPsingleKpt(FiretaskBase):
    """
    Execute IRVSP in current directory.

    """
    required_params = ["set_spn", "symprec"]
    def run_task(self, fw_spec):

        wd = os.getcwd()
        set_spn = self["set_spn"]
        symprec = self["symprec"]
        irvsp_caller = IRVSPCaller(wd, set_spn=set_spn, symprec=symprec)

        try:
            raw_struct = Structure.from_file(wd + "/POSCAR")
            formula = raw_struct.composition.formula
            structure = raw_struct.as_dict()

            outcar = Outcar(wd + "/OUTCAR")
            efermi = outcar.efermi

        except:
            formula = None
            structure = None
            efermi = None

        general = IRVSPOutputAll(wd + "/outir.txt")
        data = general.as_dict().copy()
        data["parity_eigenvals"] = {"single_kpt": general.parity_eigenvals}

        return FWAction(
            update_spec={
                "irvsp_out": data,
                "structure": structure,
                "formula": formula,
                "efermi": efermi,
                "post_relax_sg_name": irvsp_caller.sg_name,
                "post_relax_sg_number": irvsp_caller.sg_number
            }
        )

@explicit_serialize
class StandardizeCell(FiretaskBase):
    """
    Standardize cell with spglib and symprec=1e-2.

    """

    def run_task(self, fw_spec):

        wd = os.getcwd()

        struct = Structure.from_file(wd + "/POSCAR")

        numbers = [site.specie.number for site in struct]
        lattice = struct.lattice.matrix
        positions = struct.frac_coords

        if "magmom" in struct.site_properties:
            magmoms = struct.site_properties["magmom"]
            cell = (lattice, positions, numbers, magmoms)
        else:
            magmoms = None
            cell = (lattice, positions, numbers)

        lat, pos, nums = standardize_cell(cell, to_primitive=False, symprec=1e-2)

        structure = Structure(lat, nums, pos)

        if magmoms is not None:
            structure.add_site_property("magmom", magmoms)

        structure.to(fmt="poscar", filename="CONTCAR")

        return FWAction(update_spec={"structure": structure})


@explicit_serialize
class IRVSPToDb(FiretaskBase):
    """
    Stores data from outir.txt that is output by irvsp.

    required_params:
        irvsp_out (IRVSPOutput): output from IRVSP calculation.
        wf_uuid (str): unique wf id

    optional_params:
        db_file (str): path to the db file
        additional_fields (dict): dict of additional fields to add

    """

    required_params = ["irvsp_out"]
    optional_params = ["db_file", "additional_fields", "collection_name", "fw_spec_field", "task_fields_to_push"]

    def run_task(self, fw_spec):
        irvsp = self.get("irvsp_out") or fw_spec["irvsp_out"]

        irvsp = jsanitize(irvsp)

        additional_fields = self.get("additional_fields", {})
        d = additional_fields.copy()
        d["formula"] = fw_spec["formula"]
        d["efermi"] = fw_spec["efermi"]
        d["structure"] = fw_spec["structure"]
        d["irvsp"] = irvsp
        d["dir_name"] = os.getcwd()
        d["post_relax_sg_name"] = fw_spec["post_relax_sg_name"],
        d["post_relax_sg_number"] = fw_spec["post_relax_sg_number"]

        # Check for additional keys to set based on the fw_spec
        if self.get("fw_spec_field") and isinstance(self.get("fw_spec_field"), list):
            for key in self.get("fw_spec_field"):
                d.update({key: fw_spec[key]})
        # Automatically add prev fws information
        for prev_info_key in ["prev_fw_taskid", "prev_fw_db", "prev_fw_collection"]:
            if prev_info_key in fw_spec:
                d.update({prev_info_key: fw_spec[prev_info_key]})

        # store the results
        db_file = env_chk(self.get("db_file", ">>db_file<<"), fw_spec)
        if not db_file:
            with open("irvsp.json", "w") as f:
                f.write(json.dumps(d, default=DATETIME_HANDLER, indent=4))
        else:
            db = VaspCalcDb.from_db_file(db_file, admin=True)
            db.collection = db.db[self.get("collection_name", db.collection.name)]
            d.update({"db": db.db_name, "collection": db.collection.name})
            t_id = db.insert(d)
            logger.info("IRVSP calculation complete.")

        task_fields_to_push = self.get("task_fields_to_push", {}) or {}
        # pass entry information
        task_fields_to_push.update(
            {
                "prev_fw_taskid": "task_id",
                "prev_fw_db": "db",
                "prev_fw_collection": "collection"
            }
        )
        update_spec = {}
        if task_fields_to_push:
            if isinstance(task_fields_to_push, dict):
                for key, path_in_task_doc in task_fields_to_push.items():
                    if has(d, path_in_task_doc):
                        update_spec[key] = get(d, path_in_task_doc)
                    else:
                        logger.warning("Could not find {} in task document. Unable to push to next firetask/firework".format(path_in_task_doc))
            else:
                raise RuntimeError("Inappropriate type {} for task_fields_to_push. It must be a "
                                   "dictionary of format: {key: path} where key refers to a field "
                                   "in the spec and path is a full mongo-style path to a "
                                   "field in the task document".format(type(task_fields_to_push)))

        return FWAction(stored_data={"task_id": d.get("task_id", None)}, update_spec=update_spec)

