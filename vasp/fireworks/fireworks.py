import warnings

from fireworks import Firework

from pymatgen import Structure
from pymatgen.io.vasp.sets import (
    MPRelaxSet,
    MPHSERelaxSet,
    MVLScanRelaxSet,
    MVLGWSet,
    MPSOCSet,
)

from atomate.common.firetasks.glue_tasks import (
    PassCalcLocs,
    CopyFiles,
)
from atomate.vasp.config import (
    HALF_KPOINTS_FIRST_RELAX,
    RELAX_MAX_FORCE,
)
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs
from atomate.vasp.firetasks.run_calc import (
    RunVaspCustodian,
)
from atomate.vasp.firetasks.write_inputs import (
    WriteVaspFromIOSet,
    WriteVaspFromPMGObjects,
    WriteVaspStaticFromPrev,
    WriteVaspSOCFromPrev,
    ModifyIncar,
    WriteVaspNSCFFromPrev
)
from atomate.vasp.config import VASP_CMD, DB_FILE

from ..firetasks.firetasks import *
from ..firetasks.optics import *

class JOptimizeFW(Firework):
    def __init__(
            self,
            structure,
            name="structure optimization",
            vasp_input_set=None,
            vasp_cmd=VASP_CMD,
            override_default_vasp_params=None,
            ediffg=None,
            db_file=DB_FILE,
            force_gamma=True,
            job_type="double_relaxation_run",
            max_force_threshold=RELAX_MAX_FORCE,
            auto_npar=">>auto_npar<<",
            half_kpts_first_relax=HALF_KPOINTS_FIRST_RELAX,
            parents=None,
            vasptodb_kwargs=None,
            **kwargs
    ):
        """
        Optimize the given structure.

        Args:
            structure (Structure): Input structure.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use. Defaults to MPRelaxSet() if None.
            override_default_vasp_params (dict): If this is not None, these params are passed to
                the default vasp_input_set, i.e., MPRelaxSet. This allows one to easily override
                some settings, e.g., user_incar_settings, etc.
            vasp_cmd (str): Command to run vasp.
            ediffg (float): Shortcut to set ediffg in certain jobs
            db_file (str): Path to file specifying db credentials to place output parsing.
            force_gamma (bool): Force gamma centered kpoint generation
            job_type (str): custodian job type (default "double_relaxation_run")
            max_force_threshold (float): max force on a site allowed at end; otherwise, reject job
            auto_npar (bool or str): whether to set auto_npar. defaults to env_chk: ">>auto_npar<<"
            half_kpts_first_relax (bool): whether to use half the kpoints for the first relaxation
            parents ([Firework]): Parents of this particular Firework.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or MPRelaxSet(
            structure, force_gamma=force_gamma, **override_default_vasp_params
        )

        if vasp_input_set.incar["ISIF"] in (0, 1, 2, 7) and job_type == "double_relaxation":
            warnings.warn(
                "A double relaxation run might not be appropriate with ISIF {}".format(
                    vasp_input_set.incar["ISIF"]))

        vasptodb_kwargs = vasptodb_kwargs or {}
        if "additional_fields" not in vasptodb_kwargs:
            vasptodb_kwargs["additional_fields"] = {}
        vasptodb_kwargs["additional_fields"]["task_label"] = name

        t = []
        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        t.append(
            RunVaspCustodian(
                vasp_cmd=vasp_cmd,
                job_type=job_type,
                max_force_threshold=max_force_threshold,
                ediffg=ediffg,
                auto_npar=auto_npar,
                half_kpts_first_relax=half_kpts_first_relax,
            )
        )
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, **vasptodb_kwargs))
        super(JOptimizeFW, self).__init__(
            t,
            parents=parents,
            name="{}-{}".format(structure.composition.reduced_formula, name),
            **kwargs
        )

class JSelectiveOptFW(Firework):
    """
    Copy from OptimizeFW completely except adding firetask SelectiveDynmaicPoscar
    """
    def __init__(
            self,
            structure,
            name="structure optimization",
            vasp_input_set=None,
            vasp_cmd=VASP_CMD,
            override_default_vasp_params=None,
            ediffg=None,
            db_file=DB_FILE,
            force_gamma=True,
            job_type="double_relaxation_run",
            max_force_threshold=RELAX_MAX_FORCE,
            auto_npar=">>auto_npar<<",
            half_kpts_first_relax=HALF_KPOINTS_FIRST_RELAX,
            parents=None,
            vasptodb_kwargs=None,
            selective_dynamics=None,
            prev_calc_loc=True,
            **kwargs
    ):

        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or MPRelaxSet(
            structure, force_gamma=force_gamma, **override_default_vasp_params
        )

        if vasp_input_set.incar["ISIF"] in (0, 1, 2, 7) and job_type == "double_relaxation":
            warnings.warn(
                "A double relaxation run might not be appropriate with ISIF {}".format(
                    vasp_input_set.incar["ISIF"]))
        t = []
        if parents:
            if prev_calc_loc:
                t.append(
                    CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True)
                )
            mprelax_incar = MPRelaxSet(structure, force_gamma=force_gamma, **override_default_vasp_params).incar.as_dict()
            mprelax_incar.pop("@module")
            mprelax_incar.pop("@class")
            t.append(WriteVaspStaticFromPrev(other_params={"user_incar_settings": mprelax_incar}))

        elif structure:
            t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))

        if selective_dynamics:
            t.append(SelectiveDynmaicPoscar(selective_dynamics=selective_dynamics, nsites=len(structure.sites)))
        t.append(
            RunVaspCustodian(
                vasp_cmd=vasp_cmd,
                job_type=job_type,
                max_force_threshold=max_force_threshold,
                ediffg=ediffg,
                auto_npar=auto_npar,
                half_kpts_first_relax=half_kpts_first_relax,
            )
        )
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name}, **vasptodb_kwargs))
        super(JSelectiveOptFW, self).__init__(
            t,
            parents=parents,
            name="{}-{}".format(structure.composition.reduced_formula, name),
            **kwargs
        )



class JMVLGWFW(Firework):
    def __init__(
            self,
            structure,
            mode,
            name="GW",
            copy_chargcar=False,
            prev_incar=None,
            nbands=None,
            reciprocal_density=100,
            nbands_factor=5,
            ncores=16,

            vasp_input_set=None,
            vasp_input_set_params=None,
            vasp_cmd=VASP_CMD,
            prev_calc_loc=True,
            prev_calc_dir=None,
            db_file=DB_FILE,
            vasptodb_kwargs=None,
            parents=None,
            **kwargs
    ):
        """
        Standard static calculation Firework - either from a previous location or from a structure.

        Args:
            structure (Structure): Input structure. Note that for prev_calc_loc jobs, the structure
                is only used to set the name of the FW and any structure with the same composition
                can be used.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use (for jobs w/no parents)
                Defaults to MPStaticSet() if None.
            vasp_input_set_params (dict): Dict of vasp_input_set kwargs.
            vasp_cmd (str): Command to run vasp.
            prev_calc_loc (bool or str): If true (default), copies outputs from previous calc. If
                a str value, retrieves a previous calculation output by name. If False/None, will create
                new static calculation using the provided structure.
            prev_calc_dir (str): Path to a previous calculation to copy from
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            vasptodb_kwargs (dict): kwargs to pass to VaspToDb
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []

        vasp_input_set_params = vasp_input_set_params or {}
        vasptodb_kwargs = vasptodb_kwargs or {}
        if "additional_fields" not in vasptodb_kwargs:
            vasptodb_kwargs["additional_fields"] = {}
        vasptodb_kwargs["additional_fields"]["task_label"] = name

        fw_name = "{}-{}".format(
            structure.composition.reduced_formula if structure else "unknown", name
        )

        additional_file = []

        if mode == "DIAG":
            additional_file.append("WAVECAR")
            if copy_chargcar:
                additional_file.append("CHGCAR")
        elif mode == "GW":
            additional_file.append("WAVECAR")
            additional_file.append("WAVEDER")
        elif mode == "BSE":
            additional_file.append("WAVECAR")
            additional_file.append("WAVEDER")
            additional_file.append("WFULL")


        if prev_calc_dir:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_dir, contcar_to_poscar=True, additional_files=additional_file))
            t.append(WriteMVLGWFromPrev(nbands=nbands, reciprocal_density=reciprocal_density,
                                         nbands_factor=nbands_factor, ncores=ncores, prev_incar=prev_incar,
                                         mode=mode, other_params=vasp_input_set_params))
        elif parents:
            if prev_calc_loc:
                t.append(
                    CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True, additional_files=additional_file)
                )
            t.append(WriteMVLGWFromPrev(nbands=nbands, reciprocal_density=reciprocal_density,
                                         nbands_factor=nbands_factor, ncores=ncores, prev_incar=prev_incar,
                                         mode=mode, other_params=vasp_input_set_params))
        elif structure:
            vasp_input_set = vasp_input_set or MVLGWSet(
                structure, **vasp_input_set_params
            )
            t.append(
                WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
            )
        else:
            raise ValueError("Must specify structure or previous calculation")

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<", handler_group="no_handler"))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, defuse_unsuccessful="fizzle", **vasptodb_kwargs))
        super(JMVLGWFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class ScanRelaxFW(Firework):
    def __init__(
            self,
            structure,
            name="JSCAN structure optimization",
            vasp_input_set=None,
            vasp_cmd=VASP_CMD,
            override_default_vasp_params=None,
            ediffg=None,
            db_file=DB_FILE,
            force_gamma=True,
            job_type="double_relaxation_run",
            max_force_threshold=RELAX_MAX_FORCE,
            auto_npar=">>auto_npar<<",
            half_kpts_first_relax=HALF_KPOINTS_FIRST_RELAX,
            parents=None,
            vasptodb_kwargs=None,
            **kwargs
    ):
        """
        Optimize the given structure.

        Args:
            structure (Structure): Input structure.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use. Defaults to MPRelaxSet() if None.
            override_default_vasp_params (dict): If this is not None, these params are passed to
                the default vasp_input_set, i.e., MPRelaxSet. This allows one to easily override
                some settings, e.g., user_incar_settings, etc.
            vasp_cmd (str): Command to run vasp.
            ediffg (float): Shortcut to set ediffg in certain jobs
            db_file (str): Path to file specifying db credentials to place output parsing.
            force_gamma (bool): Force gamma centered kpoint generation
            job_type (str): custodian job type (default "double_relaxation_run")
            max_force_threshold (float): max force on a site allowed at end; otherwise, reject job
            auto_npar (bool or str): whether to set auto_npar. defaults to env_chk: ">>auto_npar<<"
            half_kpts_first_relax (bool): whether to use half the kpoints for the first relaxation
            parents ([Firework]): Parents of this particular Firework.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or MVLScanRelaxSet(
            structure, force_gamma=force_gamma, **override_default_vasp_params
        )

        if vasp_input_set.incar["ISIF"] in (0, 1, 2, 7) and job_type == "double_relaxation":
            warnings.warn(
                "A double relaxation run might not be appropriate with ISIF {}".format(
                    vasp_input_set.incar["ISIF"]))

        vasptodb_kwargs = vasptodb_kwargs or {}
        if "additional_fields" not in vasptodb_kwargs:
            vasptodb_kwargs["additional_fields"] = {}
        vasptodb_kwargs["additional_fields"]["task_label"] = name

        fw_name = "{}-{}".format(structure.composition.reduced_formula if structure else "unknown", name)

        t = []
        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))

        magmom = MPRelaxSet(structure).incar.get("MAGMOM", None)
        if magmom:
            t.append(ModifyIncar(incar_update={"MAGMOM": magmom}))

        t.append(
            RunVaspCustodian(
                vasp_cmd=vasp_cmd,
                job_type=job_type,
                max_force_threshold=max_force_threshold,
                ediffg=ediffg,
                auto_npar=auto_npar,
                half_kpts_first_relax=half_kpts_first_relax,
            )
        )

        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, **vasptodb_kwargs))
        super(ScanRelaxFW, self).__init__(
            t,
            parents=parents,
            name=fw_name,
            **kwargs
        )


class ScanStaticFW(Firework):
    def __init__(
            self,
            structure=None,
            name="SCAN_scf",
            vasp_input_set=None,
            vasp_input_set_params=None,
            vasp_cmd=VASP_CMD,
            force_gamma=True,
            prev_calc_loc=True,
            prev_calc_dir=None,
            db_file=DB_FILE,
            vasptodb_kwargs=None,
            parents=None,
            **kwargs
    ):
        """
        Standard static calculation Firework - either from a previous location or from a structure.

        Args:
            structure (Structure): Input structure. Note that for prev_calc_loc jobs, the structure
                is only used to set the name of the FW and any structure with the same composition
                can be used.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use (for jobs w/no parents)
                Defaults to MPStaticSet() if None.
            vasp_input_set_params (dict): Dict of vasp_input_set kwargs.
            vasp_cmd (str): Command to run vasp.
            prev_calc_loc (bool or str): If true (default), copies outputs from previous calc. If
                a str value, retrieves a previous calculation output by name. If False/None, will create
                new static calculation using the provided structure.
            prev_calc_dir (str): Path to a previous calculation to copy from
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            vasptodb_kwargs (dict): kwargs to pass to VaspToDb
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []

        vasp_input_set_params = vasp_input_set_params or {}
        vasptodb_kwargs = vasptodb_kwargs or {}
        if "additional_fields" not in vasptodb_kwargs:
            vasptodb_kwargs["additional_fields"] = {}
        vasptodb_kwargs["additional_fields"]["task_label"] = name

        fw_name = "{}-{}".format(
            structure.composition.reduced_formula if structure else "unknown", name
        )

        if prev_calc_dir:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_dir, contcar_to_poscar=True))
            t.append(WriteScanVaspStaticFromPrev(other_params=vasp_input_set_params))
        elif parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))
            t.append(WriteScanVaspStaticFromPrev(other_params=vasp_input_set_params))
        elif structure:
            vasp_input_set = vasp_input_set or "MPScanStaticSet"
            t.append(
                WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
            )
        else:
            raise ValueError("Must specify structure or previous calculation")

        t.append(RmSelectiveDynPoscar())

        magmom = MPRelaxSet(structure).incar.get("MAGMOM", None)
        if magmom:
            t.append(ModifyIncar(incar_update={"MAGMOM": magmom}))

        t.append(ModifyIncar(incar_update={"EDIFF": 1E-5}))

        if vasp_input_set_params.get("user_incar_settings", {}):
            t.append(ModifyIncar(incar_update=vasp_input_set_params.get("user_incar_settings", {})))

        if vasp_input_set_params.get("user_kpoints_settings", {}):
            t.append(WriteVaspFromPMGObjects(kpoints=vasp_input_set_params.get("user_kpoints_settings", {})))
        else:
            t.append(WriteVaspFromPMGObjects(
                kpoints=MPRelaxSet(structure=structure, force_gamma=force_gamma).kpoints.as_dict()))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, **vasptodb_kwargs))
        super(ScanStaticFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class JHSEStaticFW(Firework):
    def __init__(self, structure=None, name="HSE_scf", vasp_input_set=None, vasp_input_set_params=None,
                 vasp_cmd=VASP_CMD, prev_calc_loc=True, prev_calc_dir=None, db_file=DB_FILE, vasptodb_kwargs=None,
                 parents=None, force_gamma=True, default_magmom=True, **kwargs):
        t = []

        vasp_input_set_params = vasp_input_set_params or {}
        vasptodb_kwargs = vasptodb_kwargs or {}
        if "additional_fields" not in vasptodb_kwargs:
            vasptodb_kwargs["additional_fields"] = {}
        vasptodb_kwargs["additional_fields"]["task_label"] = name

        fw_name = "{}-{}".format(structure.composition.reduced_formula if structure else "unknown", name)


        if prev_calc_dir and parents:
            t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))
            t.append(CopyFiles(from_dir=prev_calc_dir, files_to_copy=["CHGCAR.gz", "CHGCAR", "WAVECAR.gz", "WAVECAR"],
                               continue_on_missing=True))
            t.append(WriteVaspHSEBSFromPrev(mode="uniform", reciprocal_density=None, kpoints_line_density=None))
            t.append(ModifyIncar(incar_update={"ICHARG": 11}))
        elif prev_calc_dir:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_dir, contcar_to_poscar=True, additional_files=["WAVECAR", "CHGCAR"]))
            t.append(WriteVaspHSEBSFromPrev(mode="uniform", reciprocal_density=None, kpoints_line_density=None))
        elif parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))
            t.append(WriteVaspHSEBSFromPrev(mode="uniform", reciprocal_density=None, kpoints_line_density=None))
        elif structure:
            vasp_input_set = vasp_input_set or "MPHSERelaxSet"
            # incar_hse_bs = MPHSEBSSet(structure).incar.as_dict()
            # for x in ['@module', '@class', "MAGMOM"]:
            #     incar_hse_bs.pop(x)
            t.append(
                WriteVaspFromIOSet(
                    structure=structure,
                    vasp_input_set=vasp_input_set,
                    vasp_input_params={"user_incar_settings": {
                        "EDIFF": 1E-5,
                        "NSW": 0,
                        "ISMEAR": 0,
                        "SIGMA": 0.05,
                        "ISYM": 3,
                        "NELMIN": 5
                    }}
                )
            )
        else:
            raise ValueError("Must specify structure or previous calculation")

        t.append(RmSelectiveDynPoscar())

        if default_magmom:
            magmom = MPRelaxSet(structure).incar.get("MAGMOM", None)
            t.append(ModifyIncar(incar_update={"MAGMOM": magmom}))

        if vasp_input_set_params.get("user_incar_settings", {}):
            t.append(ModifyIncar(incar_update=vasp_input_set_params.get("user_incar_settings", {})))

        if vasp_input_set_params.get("user_kpoints_settings", {}):
            t.append(WriteVaspFromPMGObjects(kpoints=vasp_input_set_params.get("user_kpoints_settings", {})))
        else:
            t.append(WriteVaspFromPMGObjects(
                kpoints=MPHSERelaxSet(structure=structure, force_gamma=force_gamma).kpoints.as_dict()))

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, **vasptodb_kwargs))
        super(JHSEStaticFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)

class JHSESOCFW(Firework):
    def __init__(
            self,
            read_chgcar=True,
            read_wavecar=False,
            magmom=None,
            structure=None,
            name="HSE_soc",
            saxis=(0, 0, 1),
            prev_calc_dir=None,
            vasp_cmd=">>vasp_ncl<<",
            copy_vasp_outputs=True,
            vasp_input_set_params = None,
            db_file=DB_FILE,
            parents=None,
            vasptodb_kwargs=None,
            **kwargs
    ):
        """
        Firework for spin orbit coupling calculation.
        Args:
            structure (Structure): Input structure. If copy_vasp_outputs, used only to set the
                name of the FW.
            name (str): Name for the Firework.
            prev_calc_dir (str): Path to a previous calculation to copy from
            vasp_cmd (str): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous
                run. Defaults to True.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        fw_name = "{}-{}".format(
            structure.composition.reduced_formula if structure else "unknown", name
        )

        vasp_input_set_params = vasp_input_set_params or {}
        vasp_input_set_params["user_incar_settings"] = {"ISYM": 3}

        vasptodb_kwargs = vasptodb_kwargs or {}
        if "additional_fields" not in vasptodb_kwargs:
            vasptodb_kwargs["additional_fields"] = {}
        vasptodb_kwargs["additional_fields"]["task_label"] = name

        if not magmom:
            magmom = [[0,0,mag_z] for mag_z in MPRelaxSet(structure).incar.get("MAGMOM", None)]

        copy_add_files_from_prev = []
        if read_chgcar:
            copy_add_files_from_prev.append("CHGCAR")
        if read_wavecar:
            copy_add_files_from_prev.append("WAVECAR")

        t = []
        if prev_calc_dir:
            t.append(
                CopyVaspOutputs(
                    calc_dir=prev_calc_dir,
                    additional_files=copy_add_files_from_prev,
                    contcar_to_poscar=True
                )
            )
            t.append(
                WriteVaspSOCFromPrev(prev_calc_dir=".", saxis=saxis)
            )
        elif parents and copy_vasp_outputs:
            t.append(
                CopyVaspOutputs(
                    calc_loc=True,
                    additional_files=copy_add_files_from_prev,
                    contcar_to_poscar=True
                )
            )
            t.append(
                WriteVaspSOCFromPrev(prev_calc_dir=".", saxis=saxis)
            )
        elif structure:
            try:
                structure.remove_site_property("magmom")
            except KeyError:
                pass
            structure.add_site_property("magmom",magmom)

            vasp_input_set = MPSOCSet(structure, nbands_factor=2)
            t.append(
                WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
            )
        else:
            raise ValueError("Must specify structure or previous calculation.")


        if vasp_input_set_params.get("user_incar_settings", {}):
            t.append(ModifyIncar(incar_update=vasp_input_set_params.get("user_incar_settings", {})))

        if vasp_input_set_params.get("user_kpoints_settings", {}):
            t.append(WriteVaspFromPMGObjects(kpoints=vasp_input_set_params.get("user_kpoints_settings", {})))

        t.extend(
            [
                RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"),
                PassCalcLocs(name=name),
                VaspToDb(db_file=db_file, **vasptodb_kwargs),
            ]
        )
        super(JHSESOCFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)

class JHSERelaxFW(Firework):
    def __init__(
            self,
            structure=None,
            name="HSE_relax",
            vasp_input_set_params=None,
            vasp_input_set=None,
            vasp_cmd=VASP_CMD,
            prev_calc_loc=True,
            prev_calc_dir=None,
            db_file=DB_FILE,
            vasptodb_kwargs=None,
            parents=None,
            force_gamma=True,
            job_type="double_relaxation_run",
            max_force_threshold=False,
            ediffg=None,
            auto_npar=">>auto_npar<<",
            default_magmom=True,
            half_kpts_first_relax=HALF_KPOINTS_FIRST_RELAX,
            **kwargs
    ):

        t = []
        vasp_input_set_params = vasp_input_set_params or {}
        vasptodb_kwargs = vasptodb_kwargs or {}
        if "additional_fields" not in vasptodb_kwargs:
            vasptodb_kwargs["additional_fields"] = {}
        vasptodb_kwargs["additional_fields"]["task_label"] = name

        fw_name = "{}-{}".format(structure.composition.reduced_formula if structure else "unknown", name)

        hse_relax_vis_incar = MPHSERelaxSet(structure=structure).incar

        if prev_calc_dir:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_dir, contcar_to_poscar=True))
            t.append(WriteVaspHSEBSFromPrev(mode="uniform", reciprocal_density=None, kpoints_line_density=None))
            t.append(ModifyIncar(incar_update=hse_relax_vis_incar))
        elif parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))
            t.append(WriteVaspHSEBSFromPrev(mode="uniform", reciprocal_density=None, kpoints_line_density=None))
            t.append(ModifyIncar(incar_update=hse_relax_vis_incar))
        elif structure:
            vasp_input_set = vasp_input_set or "MPHSERelaxSet"
            t.append(
                WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
            )
        else:
            raise ValueError("Must specify structure or previous calculation")

        if default_magmom:
            magmom = MPRelaxSet(structure).incar.get("MAGMOM", None)
            t.append(ModifyIncar(incar_update={"MAGMOM": magmom}))

        if vasp_input_set_params.get("user_incar_settings", {}):
            t.append(ModifyIncar(incar_update=vasp_input_set_params.get("user_incar_settings", {})))

        if vasp_input_set_params.get("user_kpoints_settings", {}):
            t.append(WriteVaspFromPMGObjects(kpoints=vasp_input_set_params.get("user_kpoints_settings", {})))
        else:
            t.append(WriteVaspFromPMGObjects(kpoints=MPHSERelaxSet(structure=structure,
                                                                   force_gamma=force_gamma).kpoints.as_dict()))

        t.append(
            RunVaspCustodian(
                vasp_cmd=vasp_cmd,
                job_type=job_type,
                max_force_threshold=max_force_threshold,
                ediffg=ediffg,
                auto_npar=auto_npar,
                half_kpts_first_relax=half_kpts_first_relax,
            )
        )

        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, **vasptodb_kwargs))
        super(JHSERelaxFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class JHSEcDFTFW(Firework):

    def __init__(self, up_occupation, down_occupation, nbands, prev_calc_dir=None, filesystem=None, port=27017,
                 structure=None, specific_structure=None,
                 name="HSE_cDFT", default_magmom=True,
                 vasp_input_set_params=None, job_type="normal", max_force_threshold=None,
                 vasp_cmd=VASP_CMD, db_file=DB_FILE, vasptodb_kwargs=None,
                 parents=None, prev_calc_loc=True, selective_dynamics=None, force_gamma=True, **kwargs):

        t = []
        vasp_input_set_params = vasp_input_set_params or {}
        vasptodb_kwargs = vasptodb_kwargs or {}
        if "additional_fields" not in vasptodb_kwargs:
            vasptodb_kwargs["additional_fields"] = {}
        vasptodb_kwargs["additional_fields"]["task_label"] = name

        fw_name = "{}-{}".format(structure.composition.reduced_formula if structure else "unknown", name)

        if prev_calc_dir:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_dir, additional_files=["WAVECAR"], contcar_to_poscar=True,
                                     filesystem=filesystem, port=port))
            t.append(WriteVaspHSEBSFromPrev(mode="uniform", reciprocal_density=None, kpoints_line_density=None))
            if specific_structure:
                t.append(WriteVaspFromPMGObjects(poscar=specific_structure))
            t.append(RmSelectiveDynPoscar())
        elif parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))

            t.append(WriteVaspHSEBSFromPrev(mode="uniform", reciprocal_density=None, kpoints_line_density=None))
            t.append(RmSelectiveDynPoscar())
        else:
            raise ValueError("Must specify previous calculation or parent")

        if default_magmom:
            magmom = MPRelaxSet(structure).incar.get("MAGMOM", None)
            t.append(ModifyIncar(incar_update={"MAGMOM": magmom}))

        if vasp_input_set_params.get("user_incar_settings", {}):
            t.append(ModifyIncar(incar_update=vasp_input_set_params.get("user_incar_settings", {})))
        if vasp_input_set_params.get("user_kpoints_settings", {}):
            t.append(WriteVaspFromPMGObjects(kpoints=vasp_input_set_params.get("user_kpoints_settings", {})))
        else:
            t.append(WriteVaspFromPMGObjects(
                kpoints=MPHSERelaxSet(structure=structure, force_gamma=force_gamma).kpoints.as_dict()))

        if selective_dynamics:
            t.append(SelectiveDynmaicPoscar(selective_dynamics=selective_dynamics, nsites=len(structure.sites)))

        vis_cdft = {
            "ISMEAR": -2,
            "NBANDS": nbands,
            "LDIAG": False,
            "LSUBROT": False,
            "ALGO": "All"
        }
        if up_occupation:
            vis_cdft.update({"FERWE": up_occupation})
        if down_occupation:
            vis_cdft.update({"FERDO": down_occupation})

        t.append(ModifyIncar(incar_update=vis_cdft))

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<",
                                  job_type=job_type, max_force_threshold=max_force_threshold))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, **vasptodb_kwargs))
        super(JHSEcDFTFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class JPBEcDFTRelaxFW(Firework):
    def __init__(self, prev_calc_dir, vis="MPRelaxSet", structure=None, read_structure_from=None, name="cDFT_PBE_relax",
                 vasp_input_set_params=None, vasp_cmd=VASP_CMD, db_file=DB_FILE, vasptodb_kwargs=None,
                 parents=None, wall_time=None, **kwargs):
        t = []

        vasp_input_set_params = vasp_input_set_params or {}
        vasptodb_kwargs = vasptodb_kwargs or {}
        if "additional_fields" not in vasptodb_kwargs:
            vasptodb_kwargs["additional_fields"] = {}
        vasptodb_kwargs["additional_fields"]["task_label"] = name

        fw_name = "{}-{}".format(structure.composition.reduced_formula if structure else "unknown", name)
        if read_structure_from:
            t.append(CopyVaspOutputs(additional_files=["WAVECAR"], calc_dir=read_structure_from, contcar_to_poscar=True))
        else:
            t.append(CopyVaspOutputs(additional_files=["WAVECAR"], calc_dir=prev_calc_dir))
            t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vis))
        magmom = MPRelaxSet(structure).incar.get("MAGMOM", None)
        if magmom:
            t.append(ModifyIncar(incar_update={"MAGMOM": magmom}))
        t.append(ModifyIncar(incar_update=vasp_input_set_params.get("user_incar_settings", {})))
        t.append(WriteVaspFromPMGObjects(kpoints=vasp_input_set_params.get("user_kpoints_settings", {})))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<", max_errors=5, wall_time=wall_time))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, **vasptodb_kwargs))
        super(JPBEcDFTRelaxFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class JPBEcDFTStaticFW(Firework):
    def __init__(self, structure, name="cDFT_PBE_scf",
                 vasp_input_set_params=None, vasp_cmd=VASP_CMD, db_file=DB_FILE, vasptodb_kwargs=None,
                 parents=None, wall_time=None, **kwargs):
        t = []

        vasp_input_set_params = vasp_input_set_params or {}
        vasptodb_kwargs = vasptodb_kwargs or {}
        if "additional_fields" not in vasptodb_kwargs:
            vasptodb_kwargs["additional_fields"] = {}
        vasptodb_kwargs["additional_fields"]["task_label"] = name

        fw_name = "{}-{}".format(structure.composition.reduced_formula if structure else "unknown", name)
        if parents:
            t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True))
        else:
            t.append(CopyVaspOutputs(additional_files=["CHGCAR"], calc_loc=True))
        magmom = MPRelaxSet(structure).incar.get("MAGMOM", None)
        if magmom:
            t.append(ModifyIncar(incar_update={"MAGMOM": magmom}))
        t.append(WriteVaspStaticFromPrev())
        t.append(ModifyIncar(incar_update=vasp_input_set_params.get("user_incar_settings", {})))
        t.append(WriteVaspFromPMGObjects(kpoints=vasp_input_set_params.get("user_kpoints_settings", {})))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<", max_errors=5, wall_time=wall_time))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, bandstructure_mode="uniform",
                          parse_dos=True, parse_eigenvalues=True, **vasptodb_kwargs))
        super(JPBEcDFTStaticFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)

class HSEBSFW(Firework):
    def __init__(
            self,
            parents=None,
            prev_calc_dir=None,
            cp_file_from_prev="CHGCAR",
            structure=None,
            mode="gap",
            name=None,
            input_set_overrides=None,
            vasp_cmd=VASP_CMD,
            db_file=DB_FILE,
            **kwargs
    ):
        """
        For getting a more accurate band gap or a full band structure with HSE - requires previous
        calculation that gives VBM/CBM info or the high-symmetry kpoints.

        Args:
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            prev_calc_dir (str): Path to a previous calculation to copy from
            structure (Structure): Input structure - used only to set the name of the FW.
            mode (string): options:
                "line" to get a full band structure along symmetry lines or
                "uniform" for uniform mesh band structure or
                "gap" to get the energy at the CBM and VBM
            name (str): Name for the Firework.
            vasp_cmd (str): Command to run vasp.
            db_file (str): Path to file specifying db credentials.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        name = name if name else "{} {}".format("hse", mode)

        fw_name = "{}-{}".format(
            structure.composition.reduced_formula if structure else "unknown", name
        )

        t = []
        if prev_calc_dir:
            t.append(
                CopyVaspOutputs(calc_dir=prev_calc_dir, additional_files=[cp_file_from_prev])
            )
        elif parents:
            t.append(CopyVaspOutputs(calc_loc=True, additional_files=[cp_file_from_prev]))
        else:
            raise ValueError("Must specify a previous calculation for HSEBSFW")

        t.append(WriteVaspHSEBSFromPrev(prev_calc_dir=".", mode=mode, **input_set_overrides))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))

        parse_dos = True if mode == "uniform" else False
        bandstructure_mode = mode if mode in ["line", "uniform"] else "line"

        t.append(
            VaspToDb(
                db_file=db_file,
                additional_fields={"task_label": name},
                parse_dos=parse_dos,
                bandstructure_mode=bandstructure_mode
            )
        )
        super(HSEBSFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class NonSCFFW(Firework):
    def __init__(
            self,
            parents=None,
            prev_calc_dir=None,
            filesystem=None,
            port=None,
            structure=None,
            name="nscf",
            mode="uniform",
            vasp_cmd=VASP_CMD,
            additional_files=None,
            db_file=DB_FILE,
            input_set_overrides=None,
            **kwargs
    ):
        """
        Standard NonSCF Calculation Firework supporting uniform and line modes.

        Args:
            structure (Structure): Input structure - used only to set the name
                of the FW.
            name (str): Name for the Firework.
            mode (str): "uniform" or "line" mode.
            vasp_cmd (str): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous
                run. Defaults to True.
            prev_calc_dir (str): Path to a previous calculation to copy from
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            input_set_overrides (dict): Arguments passed to the
                "from_prev_calc" method of the MPNonSCFSet. This parameter
                allows a user to modify the default values of the input set.
                For example, passing the key value pair
                    {'reciprocal_density': 1000}
                will override default k-point meshes for uniform calculations.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        input_set_overrides = input_set_overrides or {}

        fw_name = "{}-{} {}".format(
            structure.composition.reduced_formula if structure else "unknown",
            name,
            mode,
        )
        t = []

        copy_files = ["CHGCAR"]
        if additional_files:
            for file in additional_files:
                copy_files.append(file)

        if prev_calc_dir:
            t.append(
                CopyVaspOutputs(calc_dir=prev_calc_dir, additional_files=copy_files, filesystem=filesystem,
                                port=port,  contcar_to_poscar=True)
            )
        elif parents:
            t.append(CopyVaspOutputs(calc_loc=True, additional_files=copy_files))
        else:
            raise ValueError("Must specify previous calculation for NonSCFFW")

        mode = mode.lower()
        if mode == "uniform":
            t.append(
                WriteVaspNSCFFromPrev(
                    prev_calc_dir=".", mode="uniform", **input_set_overrides
                )
            )
        else:
            t.append(
                WriteVaspNSCFFromPrev(
                    prev_calc_dir=".", mode="line", **input_set_overrides
                )
            )

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(
            VaspToDb(
                db_file=db_file,
                additional_fields={"task_label": name + " " + mode},
                parse_dos=(mode == "uniform"),
                bandstructure_mode=mode,
            )
        )

        super(NonSCFFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class JSOCFW(Firework):
    def __init__(
        self,
        magmom,
        structure=None,
        name="spin-orbit coupling",
        saxis=(0, 0, 1),
        prev_calc_dir=None,
        vasp_cmd=">>vasp_ncl<<",
        copy_vasp_outputs=True,
        db_file=DB_FILE,
        parents=None,
        vasp_input_set_params=None,
        **kwargs
    ):
        """
        Firework for spin orbit coupling calculation.
        Args:
            structure (Structure): Input structure. If copy_vasp_outputs, used only to set the
                name of the FW.
            name (str): Name for the Firework.
            prev_calc_dir (str): Path to a previous calculation to copy from
            vasp_cmd (str): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous
                run. Defaults to True.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        fw_name = "{}-{}".format(
            structure.composition.reduced_formula if structure else "unknown", name
        )

        vasp_input_set_params = vasp_input_set_params or {}

        t = []
        if prev_calc_dir:
            t.append(
                CopyVaspOutputs(
                    calc_dir=prev_calc_dir,
                    additional_files=["CHGCAR"],
                    contcar_to_poscar=True,
                )
            )
            t.append(
                WriteVaspSOCFromPrev(prev_calc_dir=".", magmom=magmom, saxis=saxis, **vasp_input_set_params)
            )
        elif parents and copy_vasp_outputs:
            t.append(
                CopyVaspOutputs(
                    calc_loc=True, additional_files=["CHGCAR"], contcar_to_poscar=True
                )
            )
            t.append(
                WriteVaspSOCFromPrev(prev_calc_dir=".", magmom=magmom, saxis=saxis, **vasp_input_set_params)
            )
        elif structure:
            vasp_input_set = MPSOCSet(structure)
            t.append(
                WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
            )
        else:
            raise ValueError("Must specify structure or previous calculation.")

        t.extend(
            [
                RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"),
                PassCalcLocs(name=name),
                VaspToDb(db_file=db_file, additional_fields={"task_label": name}),
            ]
        )
        super(JSOCFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)