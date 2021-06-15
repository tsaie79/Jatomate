def rerun_irvsp(fw_id):
    from my_atomate.fireworks.pytopomat import IrvspFW
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