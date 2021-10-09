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
class RmSelectiveDynPoscar(FiretaskBase):
    def run_task(self, fw_spec):
        input_strucutre = Structure.from_file("POSCAR")
        if "selective_dynamics" in input_strucutre.site_properties.keys():
            input_strucutre.remove_site_property("selective_dynamics")
            input_strucutre.to("POSCAR", "POSCAR")

@explicit_serialize
class SelectiveDynmaicPoscar(FiretaskBase):

    required_params = ["selective_dynamics", "nsites"]

    def run_task(self, fw_spec):
        where = []
        for i in range(self["nsites"]):
            if i in self["selective_dynamics"]:
                where.append([True, True, True])
            else:
                where.append([False, False, False])
        poscar_selective = Poscar.from_file("POSCAR")
        poscar_selective.selective_dynamics = where
        poscar_selective.write_file("POSCAR")


@explicit_serialize
class WriteScanVaspStaticFromPrev(FiretaskBase):
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
        "reciprocal_density",
        "small_gap_multiply",
        "standardize",
        "sym_prec",
        "international_monoclinic",
        "lepsilon",
        "other_params",
        "potcar_spec",
    ]

    def run_task(self, fw_spec):
        lepsilon = self.get("lepsilon")

        # more k-points for dielectric calc.
        default_reciprocal_density = 200 if lepsilon else 100
        other_params = self.get("other_params", {})
        user_incar_settings = other_params.get("user_incar_settings", {})

        # for lepsilon runs, set EDIFF to 1E-5 unless user says otherwise
        if (
                lepsilon
                and "EDIFF" not in user_incar_settings
                and "EDIFF_PER_ATOM" not in user_incar_settings
        ):
            if "user_incar_settings" not in other_params:
                other_params["user_incar_settings"] = {}

        updates = {
            "ADDGRID": True,
            "LASPH": True,
            "LDAU": False,
            "LMIXTAU": True,
            "METAGGA": "SCAN",
            "NELM": 200,
        }
        other_params["user_incar_settings"].update(updates)

        vis = MPStaticSet.from_prev_calc(
            prev_calc_dir=self.get("prev_calc_dir", "."),
            reciprocal_density=self.get(
                "reciprocal_density", default_reciprocal_density
            ),
            small_gap_multiply=self.get("small_gap_multiply", None),
            standardize=self.get("standardize", False),
            sym_prec=self.get("sym_prec", 0.1),
            international_monoclinic=self.get(
                "international_monoclinic", True
            ),
            lepsilon=lepsilon,
            **other_params
        )

        potcar_spec = self.get("potcar_spec", False)
        vis.write_input(".", potcar_spec=potcar_spec)


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
        user_incar_settings = other_params.get("user_incar_settings", {})

        if "user_incar_settings" not in other_params:
            other_params["user_incar_settings"] = {}

        # updates = {
        #     # "ADDGRID": True,
        #     # "LASPH": True,
        #     # "LDAU": False,
        #     # "LMIXTAU": True,
        #     # "METAGGA": "SCAN",
        #     # "NELM": 200,
        # }
        # other_params["user_incar_settings"].update(updates)
        print(self.get("nbands"), self.get("nbands_factor"), self.get("ncores"))
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

@explicit_serialize
class FileTransferTask(FiretaskBase):
    """
    A Firetask to Transfer files.

    Before using, cp login/.ssh/id_rsa.pub to local/.ssh/authorized_keys
    then, it must already have successful scp from login to local computer, i.e.
    in OWLS: scp -P 12346 any_file jengyuantsai@localhost:any_path

    Required params:
        - mode: (str) - move, mv, copy, cp, copy2, copytree, copyfile, rtransfer
        - files: (["all"]), ([str]) or ([(str, str)]) - list of source files, or dictionary containing
                'src' and 'dest' keys
        - dest: (str) destination directory, if not specified within files parameter (else optional)

    Optional params:
        - server: (str) server host for remote transfer
        - user: (str) user to authenticate with on remote server
        - key_filename: (str) optional SSH key location for remote transfer
        - max_retry: (int) number of times to retry failed transfers; defaults to `0` (no retries)
        - retry_delay: (int) number of seconds to wait between retries; defaults to `10`
    """
    required_params = ["mode", "files", "dest", "port"]
    optional_params = ["server", "user", "key_filename", "max_retry", "retry_delay"]

    fn_list = {
        "move": shutil.move,
        "mv": shutil.move,
        "copy": shutil.copy,
        "cp": shutil.copy,
        "copy2": shutil.copy2,
        "copytree": shutil.copytree,
        "copyfile": shutil.copyfile,
    }

    def run_task(self, fw_spec):
        shell_interpret = self.get('shell_interpret', True)
        ignore_errors = self.get('ignore_errors')
        max_retry = self.get('max_retry', 0)
        retry_delay = self.get('retry_delay', 10)
        mode = self.get('mode', 'move')
        key_filename = env_chk(self.get('key_filename'), fw_spec)

        if mode == 'rtransfer':
            # remote transfers
            # Create SFTP connection
            import paramiko
            ssh = paramiko.SSHClient()
            # ssh.load_host_keys(os.path.expanduser(os.path.join("~", ".ssh", "known_hosts")))
            ssh.load_system_host_keys()
            ssh.connect(self['server'], username=self.get('user'),
                        key_filename=os.path.expanduser(os.path.join("~", ".ssh", "id_rsa")), port=self["port"])
            sftp = ssh.open_sftp()


        for f in self["files"]:
            try:
                if "all" == f:
                    src = os.getcwd().split("/")[-1]
                    dest = os.path.join(self["dest"], src)
                    # sftp.mkdir(dest)

                    for file in glob("*"):
                        try:
                            sftp.put(file, os.path.join(dest, file))
                        except FileNotFoundError:
                            sftp.mkdir(dest)
                            sftp.put(file, os.path.join(dest, file))
                else:
                    if 'src' in f:
                        src = os.path.abspath(os.path.expanduser(os.path.expandvars(f['src']))) if shell_interpret else f['src']
                    else:
                        src = abspath(os.path.expanduser(os.path.expandvars(f))) if shell_interpret else f

                    if mode == 'rtransfer':
                        dest = self['dest']
                        if os.path.isdir(src):
                            if not self._rexists(sftp, dest):
                                sftp.mkdir(dest)

                            for f in os.listdir(src):
                                if os.path.isfile(os.path.join(src,f)):
                                    sftp.put(os.path.join(src, f), os.path.join(dest, f))
                        else:
                            if not self._rexists(sftp, dest):
                                sftp.mkdir(dest)

                            sftp.put(src, os.path.join(dest, os.path.basename(src)))

                    else:
                        if 'dest' in f:
                            dest = os.path.abspath(os.path.expanduser(os.pathexpandvars(f['dest']))) if shell_interpret else f['dest']
                        else:
                            dest = os.path.abspath(os.path.expanduser(os.path.expandvars(self['dest']))) if shell_interpret else self['dest']
                        FileTransferTask.fn_list[mode](src, dest)

            except:
                traceback.print_exc()
                if max_retry:

                    # we want to avoid hammering either the local or remote machine
                    time.sleep(retry_delay)
                    self['max_retry'] -= 1
                    self.run_task(fw_spec)

                elif not ignore_errors:
                    raise ValueError(
                        "There was an error performing operation {} from {} "
                        "to {}".format(mode, self["files"], self["dest"]))

        if mode == 'rtransfer':
            sftp.close()
            ssh.close()

    @staticmethod
    def _rexists(sftp, path):
        """
        os.path.exists for paramiko's SCP object
        """
        try:
            sftp.stat(path)
        except IOError as e:
            if e[0] == 2:
                return False
            raise
        else:
            return True

@explicit_serialize
class FileSCPTask(FiretaskBase):
    """
    A Firetask to Transfer files.
    Simply using bash command SCP

    Before using, cp login/.ssh/id_rsa.pub to local/.ssh/authorized_keys
    then, it must already have successful scp from login to local computer, i.e.
    in OWLS: scp -P 12346 any_file jengyuantsai@localhost:any_path

    Required params:
        port: tunnel port
        user: user name of local workstation
        dest: absolute path in local workstation
    """
    required_params = ["port", "user", "dest"]

    def run_task(self, fw_spec):
        cmd = "scp -P {} -r {} {}@localhost:{}".format(
            self["port"],
            os.getcwd(),
            self["user"],
            self["dest"]
        )
        try:
            subprocess.check_output(cmd.split(" "))
        except subprocess.CalledProcessError as e:
            print(e.output)
            raise BaseException(e.output)


@explicit_serialize
class WriteInputsFromDB(FiretaskBase):
    """
    A Firetask to write files:
    Required params:
        - files_to_write: ([{filename:(str), contents:(str)}]) List of dicts with filenames
            and contents
    Optional params:
        - dest: (str) Shared path for files
    """
    required_params = ["db_file", "task_id", "write_chgcar"]
    optional_params = ["dest", "modify_incar"]

    def run_task(self, fw_spec):
        pth = self.get("dest", os.getcwd())
        db = VaspCalcDb.from_db_file(self["db_file"])
        e = db.collection.find_one({"task_id": self.get("task_id")})

        poscar = Poscar.from_dict(e["orig_inputs"]["poscar"])
        poscar.write_file(os.path.join(pth, "POSCAR"))

        incar = Incar.from_dict(e["orig_inputs"]["incar"])
        incar.update(self.get("modify_incar", {}))
        incar.write_file(os.path.join(pth, "INCAR"))

        kpoints = Kpoints.from_dict(e["orig_inputs"]["kpoints"])
        kpoints.write_file(os.path.join(pth, "KPOINTS"))

        if self.get("write_chgcar"):
            chgcar = db.get_chgcar(self["task_id"])
            chgcar.write_file(os.path.join(pth, "CHGCAR"))

@explicit_serialize
class Write2dNSCFKpoints(FiretaskBase):
    required_params = ["is_hse"]
    optional_params = ["added_kpoints", "reciprocal_density", "kpoints_line_density", "mode"]
    def run_task(self, fw_spec):
                 #structure, added_kpoints=None, reciprocal_density=50, kpoints_line_density=20, mode="line")
        structure = None
        try:
            structure = Structure.from_file("POSCAR")
        except Exception:
            structure = Structure.from_file("POSCAR.gz")

        added_kpoints = self.get("added_kpoints", [])
        reciprocal_density = self.get("reciprocal_density", 50)
        kpoints_line_density = self.get("kpoints_line_density", 20)
        mode = self.get("mode", "line")

        kpts = []
        weights = []
        all_labels = []

        # for both modes, include the Uniform mesh w/standard weights
        grid = Kpoints.automatic_density_by_vol(structure, reciprocal_density).kpts
        ir_kpts = SpacegroupAnalyzer(structure, symprec=0.1).get_ir_reciprocal_mesh(
            grid[0]
        )
        if self["is_hse"]:
            for k in ir_kpts:
                if round(k[0][2], 1) != 0:
                    continue
                kpts.append(k[0])
                weights.append(int(k[1]))
                all_labels.append(None)


        # for both modes, include any user-added kpoints w/zero weight
        for k in added_kpoints:
            kpts.append(k)
            weights.append(0.0)
            all_labels.append("user-defined")

        # for line mode only, add the symmetry lines w/zero weight
        if mode.lower() == "line":
            kpath = HighSymmKpath(structure)
            frac_k_points, labels = kpath.get_kpoints(
                line_density=kpoints_line_density, coords_are_cartesian=False
            )

            two_d_kpt, two_d_kpt_label = [], []
            for kpt, klabel in zip(frac_k_points, labels):
                if round(kpt[2], 1) == 0:
                    two_d_kpt.append(kpt)
                    two_d_kpt_label.append(klabel)
            frac_k_points, labels = two_d_kpt, two_d_kpt_label

            for k, f in enumerate(frac_k_points):
                kpts.append(f)
                weights.append(0.0) if self["is_hse"] else weights.append(1.0)
                all_labels.append(labels[k])

            style = Kpoints.supported_modes.Reciprocal
            kpts = kpts
            kpts_weights = weights
            labels = all_labels
            num_kpts = len(kpts)

        elif not self["is_hse"] and mode == "uniform":
            style = Kpoints.supported_modes.Gamma
            num_kpts = 0
            kpts = grid
            kpts_weights = None
            labels = None


        comment = (
            "is_HSE={} run along symmetry lines".format(self["is_hse"])
            if mode.lower() == "line"
            else "is_HSE={} run on uniform grid".format(self["is_hse"])
        )

        Kpoints(
            comment=comment,
            style=style,
            num_kpts=num_kpts,
            kpts=kpts,
            kpts_weights=kpts_weights,
            labels=labels,
        ).write_file("KPOINTS")


@explicit_serialize
class WriteVaspHSEBSFromPrev(FiretaskBase):
    """
    Writes input files for HSE band structure run. Assumes that output files
    from a previous job can be accessed. Since HSE always re-optimizes the
    charge density (no nSCF mode), the previous job is used to get the location
    of VBM/CBM for mode="gap" (otherwise just used to get the structure /
    starting charge density).

    Optional params:
        potcar_spec (bool): Instead of writing the POTCAR, write a
            "POTCAR.spec". This is intended to allow testing of workflows
            without requiring pseudo-potentials to be installed on the system.
        (documentation for all other optional params can be found in
        MPHSEBSSet)
    """

    required_params = []
    optional_params = [
        "prev_calc_dir",
        "mode",
        "reciprocal_density",
        "kpoints_line_density",
        "potcar_spec",
        "other_params"
    ]

    def run_task(self, fw_spec):
        vis = MPHSEBSSet.from_prev_calc(
            self.get("prev_calc_dir", "."),
            mode=self.get("mode", "uniform"),
            reciprocal_density=self.get("reciprocal_density", 50),
            kpoints_line_density=self.get("kpoints_line_density", 10),
            copy_chgcar=False,
            **self.get("other_params", {})
        )
        potcar_spec = self.get("potcar_spec", False)
        vis.write_input(".", potcar_spec=potcar_spec)

@explicit_serialize
class VaspToDb(FiretaskBase):
    """
    Enter a VASP run into the database. Uses current directory unless you
    specify calc_dir or calc_loc.

    Optional params:
        calc_dir (str): path to dir (on current filesystem) that contains VASP
            output files. Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        parse_dos (bool): whether to parse the DOS and store in GridFS.
            Defaults to False.
        parse_potcar_file (bool): Whether to parse the potcar file. Defaults to
            True.
        bandstructure_mode (str): Set to "uniform" for uniform band structure.
            Set to "line" for line mode. If not set, band structure will not
            be parsed.
        additional_fields (dict): dict of additional fields to add
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        fw_spec_field (str): if set, will update the task doc with the contents
            of this key in the fw_spec.
        defuse_unsuccessful (bool): this is a three-way toggle on what to do if
            your job looks OK, but is actually unconverged (either electronic or
            ionic). True -> mark job as COMPLETED, but defuse children.
            False --> do nothing, continue with workflow as normal. "fizzle"
            --> throw an error (mark this job as FIZZLED)
        task_fields_to_push (dict): if set, will update the next Firework/Firetask
            spec using fields from the task document.
            Format: {key : path} -> fw.spec[key] = task_doc[path]
            The path is a full mongo-style path so subdocuments can be referneced
            using dot notation and array keys can be referenced using the index.
            E.g "calcs_reversed.0.output.outar.run_stats"
    """
    optional_params = ["calc_dir", "calc_loc", "parse_dos", "bandstructure_mode",
                       "additional_fields", "db_file", "fw_spec_field", "defuse_unsuccessful",
                       "task_fields_to_push", "parse_chgcar", "parse_aeccar",
                       "parse_potcar_file",
                       "store_volumetric_data", "parse_eigenvalues"]

    def run_task(self, fw_spec):
        # get the directory that contains the VASP dir to parse
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])["path"]

        # parse the VASP directory
        logger.info("PARSING DIRECTORY: {}".format(calc_dir))

        drone = VaspDrone(additional_fields=self.get("additional_fields"),
                          parse_dos=self.get("parse_dos", "auto"), # JCustom
                          parse_potcar_file=self.get("parse_potcar_file", True),
                          bandstructure_mode=self.get("bandstructure_mode", False),
                          parse_chgcar=self.get("parse_chgcar", False),  # deprecated
                          parse_aeccar=self.get("parse_aeccar", False),  # deprecated
                          parse_eigenvalues=self.get("parse_eigenvalues", "auto"), # Jcustom
                          store_volumetric_data=self.get("store_volumetric_data", STORE_VOLUMETRIC_DATA))

        # assimilate (i.e., parse)
        task_doc = drone.assimilate(calc_dir)

        # Check for additional keys to set based on the fw_spec
        if self.get("fw_spec_field") and isinstance(self.get("fw_spec_field"), list):
            for key in self.get("fw_spec_field"):
                task_doc.update({key: fw_spec[key]})
        # Automatically add prev fws information
        for prev_info_key in ["prev_fw_taskid", "prev_fw_db", "prev_fw_collection"]:
            if prev_info_key in fw_spec:
                task_doc.update({prev_info_key: fw_spec[prev_info_key]})

        # get the database connection
        db_file = env_chk(self.get('db_file'), fw_spec)

        # db insertion or taskdoc dump
        if not db_file:
            with open("task.json", "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
        else:
            mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
            # Add current entry information
            task_doc.update({"db": mmdb.db_name, "collection": mmdb.collection.name})
            t_id = mmdb.insert_task(
                task_doc, use_gridfs=bool(self.get("parse_dos", False))
                                     or bool(self.get("bandstructure_mode", False))
                                     or bool(self.get("parse_eigenvalues", False))
                                     or self.get("parse_chgcar", False)  # deprecated
                                     or self.get("parse_aeccar", False)  # deprecated
                                     or bool(self.get("store_volumetric_data", STORE_VOLUMETRIC_DATA)))
            logger.info("Finished parsing with task_id: {}".format(t_id))

        defuse_children = False
        if task_doc["state"] != "successful":
            defuse_unsuccessful = self.get("defuse_unsuccessful",
                                           DEFUSE_UNSUCCESSFUL)
            if defuse_unsuccessful is True:
                defuse_children = True
            elif defuse_unsuccessful is False:
                pass
            elif defuse_unsuccessful == "fizzle":
                raise RuntimeError(
                    "VaspToDb indicates that job is not successful "
                    "(perhaps your job did not converge within the "
                    "limit of electronic/ionic iterations)!")
            else:
                raise RuntimeError("Unknown option for defuse_unsuccessful: "
                                   "{}".format(defuse_unsuccessful))

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
                    if has(task_doc, path_in_task_doc):
                        update_spec[key] = get(task_doc, path_in_task_doc)
                    else:
                        logger.warning("Could not find {} in task document. Unable to push to next firetask/firework".format(path_in_task_doc))
            else:
                raise RuntimeError("Inappropriate type {} for task_fields_to_push. It must be a "
                                   "dictionary of format: {key: path} where key refers to a field "
                                   "in the spec and path is a full mongo-style path to a "
                                   "field in the task document".format(type(task_fields_to_push)))

        return FWAction(stored_data={"task_id": task_doc.get("task_id", None)},
                        defuse_children=defuse_children, update_spec=update_spec)

