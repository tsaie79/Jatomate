"""
FWs for wflows.

"""

import warnings
import os

from fireworks import Firework

from pymatgen import Structure

from atomate.vasp.config import VASP_CMD, DB_FILE
from atomate.common.firetasks.glue_tasks import PassCalcLocs, CopyFiles
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs


from pytopomat.workflows.firetasks import (
    RunIRVSP,
    RunIRVSPAll,
    IRVSPToDb
)


class IrvspFW(Firework):
    def __init__(
            self,
            parents=None,
            structure=None,
            name="irvsp",
            run_all_kpoints=True,
            wf_uuid=None,
            db_file=DB_FILE,
            prev_calc_dir=None,
            irvsp_out=None,
            irvsptodb_kwargs=None,
            **kwargs
    ):
        """
        Run IRVSP and parse the output data. Assumes you have a previous FW with the
        calc_locs passed into the current FW.

        Args:
            structure (Structure): - only used for setting name of FW
            name (str): name of this FW
            wf_uuid (str): unique wf id
            db_file (str): path to the db file
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            prev_calc_dir (str): Path to a previous calculation to copy from
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.

        """

        fw_name = "{}-{}".format(
            structure.composition.reduced_formula if structure else "unknown", name
        )

        irvsptodb_kwargs = irvsptodb_kwargs or {}
        if "additional_fields" not in irvsptodb_kwargs:
            irvsptodb_kwargs["additional_fields"] = {}
        irvsptodb_kwargs["additional_fields"]["task_label"] = name

        t = []

        if prev_calc_dir:
            t.append(
                CopyVaspOutputs(
                    calc_dir=prev_calc_dir,
                    additional_files=["CHGCAR", "WAVECAR"],
                    contcar_to_poscar=True,
                )
            )
        elif parents:
            t.append(
                CopyVaspOutputs(
                    calc_loc=True,
                    additional_files=["CHGCAR", "WAVECAR"],
                    contcar_to_poscar=True,
                )
            )
        else:
            raise ValueError("Must specify structure or previous calculation")

        if run_all_kpoints:
            t.append(RunIRVSPAll())
        else:
            t.append(RunIRVSP())

        t.extend(
            [
                PassCalcLocs(name=name),
                IRVSPToDb(db_file=db_file, wf_uuid=wf_uuid,
                          irvsp_out=irvsp_out, **irvsptodb_kwargs),
            ]
        )

        super(IrvspFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class StandardizeFW(Firework):
    def __init__(
            self,
            parents=None,
            structure=None,
            name="standardize",
            db_file=DB_FILE,
            prev_calc_dir=None,
            vasp_cmd=None,
            **kwargs
    ):
        """
        Standardize the structure with spglib.

        Args:
            structure (Structure): pmg structure.
            name (str): name of this FW
            db_file (str): path to the db file
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            prev_calc_dir (str): Path to a previous calculation to copy from
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        fw_name = "{}-{}".format(
            structure.composition.reduced_formula if structure else "unknown", name
        )

        t = []

        if prev_calc_dir:
            t.append(
                CopyVaspOutputs(
                    calc_dir=prev_calc_dir,
                    contcar_to_poscar=True,
                )
            )
        elif parents:
            t.append(
                CopyVaspOutputs(
                    calc_loc=True,
                    contcar_to_poscar=True,
                )
            )
        else:
            raise ValueError("Must specify structure or previous calculation")

        t.extend(
            [
                StandardizeCell(),
                PassCalcLocs(name=name),
            ]
        )

        super(StandardizeFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)