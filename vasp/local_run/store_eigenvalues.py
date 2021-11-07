import json, os
from monty.json import MontyEncoder, MontyDecoder
from atomate.vasp.database import VaspCalcDb
from atomate.vasp.drones import VaspDrone

local_scf_path = "/home/qimin/tsai/sdb_tsai/Research/projects/Scan2dDefect/calc_data/scf"

mmdb = VaspCalcDb.from_db_file('/mnt/sdb/tsai/scripts/update_eigen/db.json', admin=True)

for e in mmdb.collection.find({"task_label": "SCAN_scf"})[:1]:
    t_id = e["task_id"]
    print(t_id)
    if t_id != 3257:
        continue

    path = e["dir_name"].split("/")[-1]
    drone = VaspDrone(parse_eigenvalues=True, parse_dos=True)
    task_doc = drone.assimilate(os.path.join(local_scf_path, path))

    eigenvals = {}
    dos = None
    if "calcs_reversed" in task_doc:
        for eigenvalue in ("eigenvalues", "projected_eigenvalues"):
            if eigenvalue in task_doc["calcs_reversed"][0]["output"]:  # only store idx=0 data
                eigenvals[eigenvalue] = json.dumps(task_doc["calcs_reversed"][0]["output"][eigenvalue],
                                                   cls=MontyEncoder)
                del task_doc["calcs_reversed"][0]["output"][eigenvalue]

            if "dos" in task_doc["calcs_reversed"][0]:  # only store idx=0 (last step)
                dos = json.dumps(task_doc["calcs_reversed"][0]["dos"], cls=MontyEncoder)
                del task_doc["calcs_reversed"][0]["dos"]

    if eigenvals:
        for name, data in eigenvals.items():
            data_gfs_id, compression_type = mmdb.insert_gridfs(data, "{}_fs".format(name), task_id=t_id)
            mmdb.collection.update_one(
                {"task_id": t_id}, {"$set": {"calcs_reversed.0.{}_compression".format(name): compression_type}})
            mmdb.collection.update_one({"task_id": t_id}, {"$set": {"calcs_reversed.0.{}_fs_id".format(name): data_gfs_id}})
    if dos:
        dos_gfs_id, compression_type = mmdb.insert_gridfs(
            dos, "dos_fs", task_id=t_id
        )
        mmdb.collection.update_one(
            {"task_id": t_id},
            {"$set": {"calcs_reversed.0.dos_compression": compression_type}},
        )
        mmdb.collection.update_one(
            {"task_id": t_id}, {"$set": {"calcs_reversed.0.dos_fs_id": dos_gfs_id}}
        )


