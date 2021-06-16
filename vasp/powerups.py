from vasp.firetasks.firetasks import WriteTwoDBSKpoints

from atomate.utils.utils import get_fws_and_tasks
from atomate.vasp.config import (
    VDW_KERNEL_DIR
)

from atomate.vasp.firetasks.jcustom import JFileTransferTask, JWriteInputsFromDB
from atomate.vasp.firetasks.glue_tasks import CopyFiles
from atomate.vasp.firetasks.write_inputs import ModifyIncar, WriteVaspFromPMGObjects

from pymatgen import Structure
from pymatgen.io.vasp.sets import MPRelaxSet

__author__ = "Jeng-Yuan Tsai"
__email__ = "tsaie79@gmail.com"


def scp_files(
        original_wf,
        dest,
        fw_name_constraint=None,
        task_name_constraint="VaspToDb",
):
    """
    SCP ALL files to local computer

    Args:
        original_wf (Workflow)
        dest (str): "/home/jengyuantsai/test_scp_fw/defect_db/binary_vac_AB/" (make sure every folder exists)
        fw_name_constraint (str): pattern for fireworks to clean up files after
        task_name_constraint (str): pattern for firetask to clean up files

    Returns:
       Workflow
    """
    idx_list = get_fws_and_tasks(
        original_wf,
        fw_name_constraint=fw_name_constraint,
        task_name_constraint=task_name_constraint,
    )
    for idx_fw, idx_t in idx_list:
        original_wf.fws[idx_fw].tasks.insert(idx_t + 1, JFileTransferTask(
            mode="rtransfer",
            files=["all"],
            dest=dest,
            server="localhost",
            user="jengyuantsai"
        ))

    return original_wf

def write_inputs_from_db(original_wf, db_file, task_id, modify_incar, write_chgcar=True, fw_name_constraint=None):

    idx_list = get_fws_and_tasks(
        original_wf,
        fw_name_constraint=fw_name_constraint,
        task_name_constraint="RunVasp",
    )
    for idx_fw, idx_t in idx_list:
        original_wf.fws[idx_fw].tasks.insert(idx_t - 1, JWriteInputsFromDB(db_file=db_file, task_id=task_id,
                                                                           write_chgcar=write_chgcar,
                                                                           modify_incar=modify_incar))
    return original_wf

def jmodify_to_soc(
        original_wf,
        structure,
        nbands=None,
        saxis=[0,0,1],
        magmom=None,
        modify_incar_params=None,
        fw_name_constraint=None,
):
    """
    Takes a regular workflow and transforms its VASP fireworkers that are
    specified with fw_name_constraints to non-collinear calculations taking spin
    orbit coupling into account.
    Args:
        original_wf (Workflow): The original workflow.
        nbands (int): number of bands selected by the user (for now)
        structure (Structure)
        modify_incar_params ({}): a dictionary containing the setting for
            modifying the INCAR (e.g. {"ICHARG": 11})
        fw_name_constraint (string): name of the fireworks to be modified (all
            if None is passed)
    Returns:
        Workflow: modified with SOC
    """

    if structure is None:
        try:
            sid = get_fws_and_tasks(
                original_wf,
                fw_name_constraint="structure optimization",
                task_name_constraint="WriteVasp",
            )
            fw_id = sid[0][0]
            task_id = sid[0][1]
            structure = (
                original_wf.fws[fw_id].tasks[task_id]["vasp_input_set"].structure
            )
        except:
            raise ValueError(
                "modify_to_soc powerup requires the structure in vasp_input_set"
            )

    if not magmom:
        magmom = [[0,0,mag_z] for mag_z in MPRelaxSet(structure).incar.get("MAGMOM", None)]

    modify_incar_soc = {
        "incar_update": {
            "LSORBIT": "T",
            "SAXIS": saxis,
            "MAGMOM": magmom,
            "ISPIN": 2,
            "ICHARG":11,
            # "LMAXMIX": 4,
            "ISYM": 0,
        }
    }
    if nbands:
        modify_incar_soc["incar_update"].update({"NBANDS": nbands})

    if modify_incar_params:
        modify_incar_soc["incar_update"].update(modify_incar_params)

    run_vasp_list = get_fws_and_tasks(
        original_wf,
        fw_name_constraint=fw_name_constraint,
        task_name_constraint="RunVasp",
    )
    for idx_fw, idx_t in run_vasp_list:
        original_wf.fws[idx_fw].tasks[idx_t]["vasp_cmd"] = ">>vasp_ncl<<"
        original_wf.fws[idx_fw].tasks.insert(idx_t, ModifyIncar(**modify_incar_soc))

        original_wf.fws[idx_fw].name += "_soc"

    run_boltztrap_list = get_fws_and_tasks(
        original_wf,
        fw_name_constraint=fw_name_constraint,
        task_name_constraint="RunBoltztrap",
    )
    for idx_fw, idx_t in run_boltztrap_list:
        original_wf.fws[idx_fw].name += "_soc"

    # revise task_label to xxx_soc in db
    to_db_list = get_fws_and_tasks(
        original_wf,
        fw_name_constraint=fw_name_constraint,
        task_name_constraint="VaspToDb",
    )
    for idx_fw, idx_t in to_db_list:
        original_wf.fws[idx_fw].tasks[idx_t]["additional_fields"].update({"task_label": original_wf.fws[idx_fw].name})

    return original_wf

def remove_todb(original_wf, fw_name_constraint=None):
    """
    Simple powerup that clears the VaspToDb to a workflow.
    Args:
        original_wf (Workflow): The original workflow.
        fw_name_constraint (str): name constraint for fireworks to
            have their modification tasks removed
    """
    idx_list = get_fws_and_tasks(
        original_wf,
        fw_name_constraint=fw_name_constraint,
        task_name_constraint="VaspToDb",
    )
    idx_list.reverse()
    for idx_fw, idx_t in idx_list:
        original_wf.fws[idx_fw].tasks.pop(idx_t)
    return original_wf

def write_PMGObjects(original_wf, pmg_objs, fw_name_constraint=None):

    idx_list = get_fws_and_tasks(
        original_wf,
        fw_name_constraint=fw_name_constraint,
        task_name_constraint="RunVasp",
    )

    for idx_fw, idx_t in idx_list:
        original_wf.fws[idx_fw].tasks.insert(
            idx_t, WriteVaspFromPMGObjects(**pmg_objs)
        )
    return original_wf

def cp_vdw_file(original_wf, fw_name_constraint=None):

    idx_list = get_fws_and_tasks(
        original_wf,
        fw_name_constraint=fw_name_constraint,
        task_name_constraint="RunVasp",
    )

    for idx_fw, idx_t in idx_list:
        original_wf.fws[idx_fw].tasks.insert(
            idx_t-1, CopyFiles(from_dir=VDW_KERNEL_DIR)
        )
    return original_wf

def cp_vasp_from_prev(original_wf, vasp_io, fw_name_constraint=None):
    idx_list = get_fws_and_tasks(
        original_wf,
        fw_name_constraint=fw_name_constraint,
        task_name_constraint="CopyVaspOutputs",
    )
    for idx_fw, idx_t in idx_list:
        if original_wf.fws[idx_fw].tasks[idx_t]["additional_files"]:
            original_wf.fws[idx_fw].tasks[idx_t]["additional_files"].extend(vasp_io)
        else:
            original_wf.fws[idx_fw].tasks[idx_t].update({"additional_files": vasp_io})
    return original_wf

def add_modify_twod_bs_kpoints(
        original_wf, modify_kpoints_params=None, fw_name_constraint=None
):
    """
    Every FireWork that runs VASP has a ModifyKpoints task just beforehand. For
    example, allows you to modify the KPOINTS based on the Worker using env_chk
    or using hard-coded changes.

    Args:
        original_wf (Workflow)
        modify_kpoints_params (dict): dict of parameters for ModifyKpoints.
        fw_name_constraint (str): Only apply changes to FWs where fw_name
        contains this substring.

    Returns:
       Workflow
    """
    modify_kpoints_params = modify_kpoints_params or {
        "twod_kpoints_update": ">>twod_kpoints_update<<"
    }

    idx_list = get_fws_and_tasks(
        original_wf,
        fw_name_constraint=fw_name_constraint,
        task_name_constraint="RunVasp",
    )
    for idx_fw, idx_t in idx_list:
        original_wf.fws[idx_fw].tasks.insert(
            idx_t,
            WriteTwoDBSKpoints(**modify_kpoints_params)
        )
    return original_wf
