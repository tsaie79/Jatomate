import os.path

from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MPScanRelaxSet

from atomate.vasp.fireworks.core import HSEBSFW
from atomate.vasp.powerups import *
from atomate.vasp.workflows.base.core import get_wf

from ..fireworks.fireworks import *
import my_atomate

from fireworks import Workflow

import numpy as np

def get_wf_full_hse(structure, charge_states, gamma_only, gamma_mesh, nupdowns, task,
                    vasptodb=None, wf_addition_name=None, task_arg=None):

    encut = 1.3*max([potcar.enmax for potcar in MPHSERelaxSet(structure).potcar])

    print("SET ENCUT:{}".format(encut))

    vasptodb = vasptodb or {}
    task_arg = task_arg or {}

    fws = []
    for cs, nupdown in zip(charge_states, nupdowns):
        print("Formula: {}".format(structure.formula))
        if structure.site_properties.get("magmom", None):
            structure.remove_site_property("magmom")
        structure.set_charge(cs)
        nelect = MPHSERelaxSet(structure, use_structure_charge=True).nelect
        user_incar_settings = {
            "ENCUT": encut,
            "ISIF": 2,
            "ISMEAR": 0,
            "EDIFF": 1e-5,
            "EDIFFG": -0.02,
            "LCHARG": False,
            "NUPDOWN": nupdown,
            "SIGMA": 0.01,
            "NSW": 150
            #"NCORE": 4 owls normal 14; cori 8. Reduce ncore if want to increase speed but low memory risk
        }

        user_incar_settings.update({"NELECT": nelect})

        if gamma_only is True:
            # user_kpoints_settings = Kpoints.gamma_automatic((1,1,1), (0.333, 0.333, 0))
            user_kpoints_settings = Kpoints.gamma_automatic()

        elif gamma_only:
            nkpoints = len(gamma_only)
            kpts_weights = [1.0 for i in np.arange(nkpoints)]
            labels = [None for i in np.arange(nkpoints)]
            user_kpoints_settings = Kpoints.from_dict(
                {
                    'comment': 'JCustom',
                    'nkpoints': nkpoints,
                    'generation_style': 'Reciprocal',
                    'kpoints': gamma_only,
                    'usershift': (0, 0, 0),
                    'kpts_weights': kpts_weights,
                    'coord_type': None,
                    'labels': labels,
                    'tet_number': 0,
                    'tet_weight': 0,
                    'tet_connections': None,
                    '@module': 'pymatgen.io.vasp.inputs',
                    '@class': 'Kpoints'
                }
            )

        else:
            user_kpoints_settings = None


        # FW1 Structure optimization firework
        opt = JOptimizeFW(
            structure=structure,
            name="PBE_relax",
            max_force_threshold=False,
            job_type="double_relaxation_run",
            ediffg=-0.01,
            force_gamma=gamma_mesh,
            vasptodb_kwargs={
                "parse_dos": False,
                "parse_eigenvalues": False,
            },
            override_default_vasp_params={
                "user_incar_settings": user_incar_settings,
                "user_kpoints_settings": user_kpoints_settings
            },
        )

        # FW2 Run HSE relax
        def hse_relax(parents):
            fw = JHSERelaxFW(
                structure=structure,
                force_gamma=gamma_mesh,
                job_type="double_relaxation_run",
                ediffg=-0.01,
                vasp_input_set_params={
                    "user_incar_settings": user_incar_settings,
                    "user_kpoints_settings": user_kpoints_settings
                },
                name="HSE_relax",
                vasptodb_kwargs={
                    "additional_fields": {
                        "charge_state": cs,
                        "nupdown_set": nupdown
                    },
                    "parse_dos": False,
                    "parse_eigenvalues": False
                },
                parents=parents
            )
            return fw

        # FW3 Run HSE SCF
        uis_hse_scf = {
            "user_incar_settings": {
                "LVHAR": True,
                # "AMIX": 0.2,
                # "AMIX_MAG": 0.8,
                # "BMIX": 0.0001,
                # "BMIX_MAG": 0.0001,
                "EDIFF": 1E-05,
                "ENCUT": encut,
                "ISMEAR": 0,
                "LCHARG": False,
                "LWAVE": True,
                "NSW": 0,
                "NUPDOWN": nupdown,
                "NELM": 150,
                "SIGMA": 0.05
            },
            "user_kpoints_settings": user_kpoints_settings
        }

        uis_hse_scf["user_incar_settings"].update({"NELECT": nelect})

        def hse_scf(parents, prev_calc_dir=None, lcharg=False, parse_dos=True, parse_eigenvalues=True,
                    metadata_to_pass=None):

            bandstructure_mode = None
            if parse_dos:
                uis_hse_scf["user_incar_settings"].update({"ENMAX": 10, "ENMIN": -10, "NEDOS": 9000})
                bandstructure_mode = "uniform"
            else:
                bandstructure_mode = False

            if lcharg:
                uis_hse_scf["user_incar_settings"].update({"LCHARG":True})



            fw = JHSEStaticFW(
                structure,
                force_gamma=gamma_mesh,
                vasp_input_set_params=uis_hse_scf,
                prev_calc_dir=prev_calc_dir,
                parents=parents,
                name="HSE_scf",
                vasptodb_kwargs={
                    "additional_fields": {
                        "task_type": "JHSEStaticFW",
                        "charge_state": cs,
                        "nupdown_set": nupdown
                    },
                    "parse_dos": parse_dos,
                    "parse_eigenvalues": parse_eigenvalues,
                    "bandstructure_mode": bandstructure_mode,
                    "task_fields_to_push": metadata_to_pass
                },
            )
            return fw

        def hse_soc(parents, prev_calc_dir=None, parse_dos=True,
                    parse_eigenvalues=True, read_chgcar=True, read_wavecar=True, saxis=(0,0,1),
                    metadata_from_scf=None):

            if parse_dos:
                uis_hse_scf["user_incar_settings"].update({"ENMAX": 10, "ENMIN": -10, "NEDOS": 9000})
                bandstructure_mode = "uniform"

            fw = JHSESOCFW(
                prev_calc_dir=prev_calc_dir,
                structure=structure,
                read_chgcar=read_chgcar,
                read_wavecar=read_wavecar,
                name="HSE_soc",
                saxis=saxis,
                parents=parents,
                vasp_input_set_params=uis_hse_scf,
                vasptodb_kwargs={
                    "additional_fields": {
                        "task_type": "JHSESOCFW",
                        "charge_state": cs,
                        "nupdown_set": nupdown
                    },
                    "parse_dos": parse_dos,
                    "parse_eigenvalues": parse_eigenvalues,
                    "bandstructure_mode": bandstructure_mode,
                    "fw_spec_field": list(metadata_from_scf.keys()) if metadata_from_scf else []
                }
            )
            return fw


        def hse_bs(parents, mode="line", prev_calc_dir=None):
            if mode == "uniform":
                uis_hse_scf["user_incar_settings"].update({"ENMAX": 10, "ENMIN": -10, "NEDOS": 9000})

            fw = HSEBSFW(
                structure=structure,
                mode=mode,
                input_set_overrides={"other_params": {"two_d_kpoints": True,
                                                      "user_incar_settings":uis_hse_scf["user_incar_settings"],
                                                      },
                                     "kpoints_line_density": 20
                                     },
                cp_file_from_prev="CHGCAR",
                prev_calc_dir=prev_calc_dir,
                parents=parents,
                name="HSE_bs"
            )
            return fw

        if task == "opt":
            fws.append(opt)
        elif task == "hse_relax":
            fws.append(hse_relax(parents=None))
        elif task == "hse_scf":
            fws.append(hse_scf(parents=None, **task_arg))
        elif task == "hse_bs":
            fws.append(hse_bs(parents=None, **task_arg))
        elif task == "hse_soc":
            fws.append(hse_soc(parents=None, **task_arg))
        elif task == "hse_scf-hse_bs":
            fws.append(hse_scf(parents=None))
            fws.append(hse_bs(parents=fws[-1], **task_arg))
        elif task == "hse_scf-hse_soc":
            fws.append(hse_scf(parents=None, lcharg=True, **task_arg))
            fws.append(hse_soc(parents=fws[-1]))
        elif task == "hse_relax-hse_scf":
            fws.append(hse_relax(parents=None))
            fws.append(hse_scf(fws[-1], **task_arg))
        elif task == "opt-hse_relax-hse_scf":
            fws.append(opt)
            fws.append(hse_relax(parents=fws[-1]))
            fws.append(hse_scf(parents=fws[-1], **task_arg))
        elif task == "hse_relax-hse_scf-hse_bs":
            fws.append(hse_relax(parents=None))
            fws.append(hse_scf(parents=fws[-1], lcharg=True))
            fws.append(hse_bs(parents=fws[-1], **task_arg))
        elif task == "hse_relax-hse_scf-hse_soc":
            fws.append(hse_relax(parents=None))
            fws.append(hse_scf(parents=fws[-1], lcharg=True, **task_arg))
            if "metadata_to_pass" in task_arg:
                metadata_from_scf = task_arg["metadata_to_pass"]
            else:
                metadata_from_scf = None
            fws.append(hse_soc(parents=fws[-1], metadata_from_scf=metadata_from_scf))
        elif task == "opt-hse_relax-hse_scf-hse_bs":
            fws.append(opt)
            fws.append(hse_relax(parents=fws[-1]))
            fws.append(hse_scf(parents=fws[-1], lcharg=True))
            fws.append(hse_bs(parents=fws[-1], **task_arg))


    wf_name = "{}:{}:q{}:sp{}".format("".join(structure.formula.split(" ")), wf_addition_name, charge_states, nupdowns)

    wf = Workflow(fws, name=wf_name)

    vasptodb.update({"wf": [fw.name for fw in wf.fws]})
    wf = add_additional_fields_to_taskdocs(wf, vasptodb)
    wf = add_namefile(wf)
    return wf


def get_wf_full_scan(structure, charge_states, gamma_only, dos, nupdowns,
                     vasptodb=None, wf_addition_name=None, wf_yaml=None):

    encut = 1.3*max([potcar.enmax for potcar in MPScanRelaxSet(structure).potcar])
    print("SET ENCUT:{}".format(encut))

    vasptodb = vasptodb or {}

    fws = []
    for cs, nupdown in zip(charge_states, nupdowns):
        print("Formula: {}".format(structure.formula))
        if structure.site_properties.get("magmom", None):
            structure.remove_site_property("magmom")
        structure.set_charge(cs)
        nelect = MPRelaxSet(structure, use_structure_charge=True).nelect

        if gamma_only is True:
            kpt = Kpoints.gamma_automatic()
            user_kpoints_settings = kpt.__dict__
            user_kpoints_settings.update({"style": "G"})
            user_kpoints_settings.pop("_style")
            # user_kpoints_settings = MPRelaxSet(structure).kpoints

        elif gamma_only:
            nkpoints = len(gamma_only)
            kpts_weights = [1.0 for i in np.arange(nkpoints)]
            labels = [None for i in np.arange(nkpoints)]
            kpt = Kpoints.from_dict(
                {
                    'comment': 'JCustom',
                    'nkpoints': nkpoints,
                    'generation_style': 'Reciprocal',
                    'kpoints': gamma_only,
                    'usershift': (0, 0, 0),
                    'kpts_weights': kpts_weights,
                    'coord_type': None,
                    'labels': labels,
                    'tet_number': 0,
                    'tet_weight': 0,
                    'tet_connections': None,
                    '@module': 'pymatgen.io.vasp.inputs',
                    '@class': 'Kpoints'
                }
            )
            user_kpoints_settings = kpt.__dict__

        else:
            user_kpoints_settings = None

        uis = {
            "user_incar_settings": {
                "ENCUT": encut,
                "NUPDOWN": nupdown,
                "NELECT": nelect
            },
            "user_kpoints_settings": user_kpoints_settings
        }

        wf_yaml = wf_yaml if wf_yaml else os.path.join(os.path.dirname(os.path.abspath(__file__)), "general/scan.yaml")
        wf = get_wf(structure, wf_yaml)
        if dos:
            wf = add_modify_incar(wf, {"incar_update": {"EMAX": 10, "EMIN": -10, "NEDOS": 9000}}, fw_name_constraint="SCAN_scf")
        if uis.get("user_incar_settings"):
            wf = add_modify_incar(wf, {"incar_update": uis["user_incar_settings"]})
        if uis.get("user_kpoints_settings"):
            wf = add_modify_kpoints(wf, {"kpoints_update": uis["user_kpoints_settings"]})

        vasptodb.update({"wf": [fw.name for fw in wf.fws], "charge_state": cs, "nupdown_set": nupdown})
        wf = add_additional_fields_to_taskdocs(wf, vasptodb)
        wf = add_additional_fields_to_taskdocs(wf, {"charge_state": cs}, task_name_constraint="IRVSPToDb")
        fws.extend(wf.fws)

    wf_name = "{}:{}:q{}:sp{}".format("".join(structure.formula.split(" ")), wf_addition_name, charge_states, nupdowns)
    wf = Workflow(fws, name=wf_name)
    wf = add_namefile(wf)
    return wf

