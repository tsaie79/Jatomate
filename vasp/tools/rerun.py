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
def rerun_fw(fworker):
    """
    fworker = "nersc", "owls", "efrc"
    """
    from fireworks import LaunchPad
    import os, subprocess, shutil, datetime
    from glob import glob
    from pymatgen.io.vasp.inputs import Structure, Incar
    from pymatgen.io.vasp.outputs import Oszicar
    from qubitPack.tool_box import get_db

    db = get_db("mgb2", "workflows", port=12345, user="Jeng")
    lpad = LaunchPad(host=db.host, port=db.port, name=db.db_name, username=db.user, password=db.password)

    def delet_dir(dir_name):
        shutil.rmtree(dir_name)
    def rerun_opt(fw_id): #421
        # for f in glob("POSCAR*"):
        #     os.remove(f)
        lpad.rerun_fw(fw_id, recover_launch="last", recover_mode="prev_dir")
        fw = lpad.get_fw_by_id(fw_id)
        # fw.spec['_queueadapter'] = {'walltime': '06:00:00'}
        # fw.tasks[1]["other_params"]["user_incar_settings"].update({"ALGO":"All"})
        # fw.tasks[1]["incar_update"].update({"ALGO": "All"})
        # lpad.update_spec([fw.fw_id], fw.as_dict()["spec"])
        shutil.copy("CONTCAR", "POSCAR")
        # subprocess.call("qlaunch -c {} -r singleshot -f {}".format(
        #     os.path.join(os.path.expanduser("~"), "config/project/Scan2dDefect/calc_data/"), fw.fw_id).split(" "))
    def rerun_scf(fw_id):
        subprocess.call("gunzip INCAR.gz".split(" "))
        incar = Incar.from_file("INCAR")
        incar.update({"NELM": 200, "ALGO": "Fast"})
        incar.write_file("INCAR")
        lpad.rerun_fw(fw_id, recover_launch="last", recover_mode="prev_dir")




    fws = lpad.get_fw_ids(
        {
            "state": {"$in": ["READY"]}, #First running, then fizzled
            "name": {"$regex": "static"},
            "spec._fworker": fworker,
            "created_on": {"$gte": str(datetime.datetime(2021,8,12))}
            # "fw_id": {"$gte": 1797}
            # "fw_id": 2550
        }
    )
    print(fws)
    a = []
    for fw_id in fws[::-1][:]:
        prev_fw = lpad.get_fw_by_id(fw_id)
        fworker = prev_fw.spec["_fworker"]
        prev_path = lpad.get_launchdir(fw_id, 0)
        os.chdir(prev_path)
        # a.append(fw_id)
        # if glob("CONTCAR.relax*") != []:
        try:
            print(fw_id, fworker, prev_path, prev_fw.name)
            # delet_dir(prev_path)
            # rerun_opt(fw_id)
            # rerun_scf(fw_id)
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

rerun_fw("efrc")

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

#%%
from fireworks import LaunchPad
import os
lpad = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/Scan2dDefect/calc_data/my_launchpad.yaml"))
for i in lpad.get_wf_ids({"state": "COMPLETED"}):
    print(i)
    try:
        fw = lpad.get_fw_by_id(i)
        if fw.spec["_fworker"] == "owls":
            print(i, fw.spec["_fworker"])
            lpad.delete_wf(i, True)
    except Exception:
            continue
#%%
from fireworks import LaunchPad
import os
lpad = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/Scan2dMat/calc_data/my_launchpad.yaml"))

for i in lpad.get_wf_ids({"state": "READY"})[:10]:
    wf = lpad.get_wf_by_fw_id(i)
    for fw in wf.fws:
        if fw.spec["_fworker"] == "owls":
            fw.spec["_fworker"] = "efrc"
            lpad.update_spec([fw.fw_id], fw.as_dict()["spec"])


#%%
from fireworks import LaunchPad
import os

lpad = LaunchPad.from_file(
    os.path.join(os.path.expanduser("~"), "config/project/Scan2dDefect/calc_data/my_launchpad.yaml"))
for i in lpad.get_wf_ids({"state": "COMPLETED"}):
    print(i)
    try:
        fw = lpad.get_fw_by_id(i)
        if fw.spec["_fworker"] == "owls":
            print(i, fw.spec["_fworker"])
            lpad.delete_wf(i, True)
    except Exception:
        continue




