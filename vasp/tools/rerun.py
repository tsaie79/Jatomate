#%%
def rerun_irvsp(fw_id):
    from vasp.fireworks import IrvspFW
    import os, subprocess

    p = lpad.get_launchdir(fw_id)
    line = lpad.get_fw_by_id(fw_id+1)
    prev = line.launches[-1].launch_dir
    c2db_uid = line.tasks[-1]["additional_fields"]["c2db_uid"]
    ff = IrvspFW(
        prev_calc_dir=prev,
        symprec=0.001,
        irvsptodb_kwargs=dict(collection_name="ir_data", additional_fields={"c2db_uid": c2db_uid})
    )
    lpad.rerun_fw(fw_id)
    lpad.update_spec([fw_id], ff.as_dict()["spec"])
    os.chdir(p)
    subprocess.call("rlaunch -c /global/homes/t/tsaie79/config/project/C2DB_IR/calc_data/ singleshot -f {}".format(fw_id).split(" "))


#%%
def rerun_relax(fworker):
    """
    fworker = "nersc", "owls", "efrc"
    """
    from fireworks import LaunchPad
    import os, subprocess, shutil
    from glob import glob
    from pymatgen.io.vasp.inputs import Structure
    from pymatgen.io.vasp.outputs import Oszicar

    def delet_dir(fw):
        dir_name = lpad.get_launchdir(fw.fw_id)
        shutil.rmtree(dir_name)
    def rerun_opt(fw): #421
        # for f in glob("POSCAR*"):
        #     os.remove(f)
        # subprocess.call("gunzip OUTCAR.gz CONTCAR.gz".split(" "))
        # st = Structure.from_file("CONTCAR")
        # fw.tasks[0].update({"structure": st})
        fw.spec["_queueadapter"].update({"walltime": "06:00:00"})
        fw.tasks[2]["incar_update"].update({"ISPIN": 1})
        lpad.rerun_fw(fw.fw_id)
        lpad.update_spec([fw.fw_id], fw.as_dict()["spec"])
        # subprocess.call("qlaunch -c {} -r singleshot -f {}".format(
        #     os.path.join(os.path.expanduser("~"), "config/project/Scan2dMat/calc_data/"), fw.fw_id).split(" "))

    lpad = LaunchPad.from_file(
        os.path.join(os.path.expanduser("~"), "config/project/Scan2dMat/calc_data/my_launchpad.yaml"))

    fws = lpad.get_fw_ids(
        {
            "state": {"$in": ["FIZZLED"]}, #First running, then fizzled
            "name": {"$regex": "SCAN_relax pos"},
            "spec._fworker": fworker,
            # "fw_id": {"$nin": [1896, 2160, 2172, 2484, 138]}
            # "fw_id": 2550
        }
    )
    a = []
    for fw_id in fws:
        fw = lpad.get_fw_by_id(fw_id)
        fworker = fw.spec["_fworker"]
        path = lpad.get_launchdir(fw_id)
        # os.chdir(path)
        # a.append(fw_id)
        # if glob("CONTCAR.relax*") != []:
        try:
            print(fw_id, fworker, path, fw.name)
            rerun_opt(fw)
            # delet_dir(fw)
        except Exception as err:
            print(err)
            continue
        # try:
        #     oszicar = Oszicar("OSZICAR")
        # except Exception as err:
        #     print(fw_id, err)
        #     continue
        # if len(Oszicar("OSZICAR").ionic_steps) > 10:
        #     a.append(fw_id)
        #     print(fw_id, fworker, path, fw.name)
        #     # rerun_opt(fw)
    print(a, len(a))

rerun_relax("owls")

#%%
def rerun_irvsp():
    lpad = LaunchPad.from_file(
        os.path.join(os.path.expanduser("~"), "config/project/Scan2dMat/calc_data/my_launchpad.yaml"))

    for fw_id in lpad.get_fw_ids(
            {"state": {"$in": ["RUNNING"]}, "name": {"$regex":"irvsp"}, "spec._fworker": fworker}
    ):
        fw = lpad.get_fw_by_id(fw_id)
        fworker = fw.spec["_fworker"]
        path = lpad.get_launchdir(fw_id)
        os.chdir(path)
        print(fw_id, fworker, path, fw.name)

        fw.spec["_queueadapter"].update({"walltime": "01:00:00"})
        lpad.rerun_fw(fw_id)
        lpad.update_spec([fw.fw_id], fw.as_dict()["spec"])
        subprocess.call("qlaunch -c {} -r singleshot -f {}".format(
            os.path.join(os.path.expanduser("~"), "config/project/Scan2dMat/calc_data/"), fw_id).split(" "))



