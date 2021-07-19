#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# This file has been created at the Institute of
# Landscape Systems Analysis at the ZALF.
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

import capnp
from collections import defaultdict
import copy
import csv
from datetime import date, timedelta
import json
import math
import numpy as np
import os
from pathlib import Path
from pyproj import CRS, Transformer
import sys
import time

import monica_io3
import soil_io3
import monica_run_lib as Mrunlib

PATH_TO_REPO = Path(os.path.realpath(__file__)).parent
if str(PATH_TO_REPO) not in sys.path:
    sys.path.insert(1, str(PATH_TO_REPO))

PATH_TO_MAS_INFRASTRUCTURE_REPO = PATH_TO_REPO.parent / "mas-infrastructure"
if str(PATH_TO_MAS_INFRASTRUCTURE_REPO) not in sys.path:
    sys.path.insert(1, str(PATH_TO_MAS_INFRASTRUCTURE_REPO))

import src.python.common.common as common

PATH_TO_CAPNP_SCHEMAS = PATH_TO_MAS_INFRASTRUCTURE_REPO / "capnproto_schemas"
abs_imports = [str(PATH_TO_CAPNP_SCHEMAS)]
soil_data_capnp = capnp.load(str(PATH_TO_CAPNP_SCHEMAS / "soil_data.capnp"), imports=abs_imports) 
common_capnp = capnp.load(str(PATH_TO_CAPNP_SCHEMAS / "common.capnp"), imports=abs_imports) 
climate_data_capnp = capnp.load(str(PATH_TO_CAPNP_SCHEMAS / "climate_data.capnp"), imports=abs_imports)
model_capnp = capnp.load(str(PATH_TO_CAPNP_SCHEMAS / "model.capnp"), imports=abs_imports)

PATHS = {
    # adjust the local path to your environment
    "mbm-local-remote": {
        "include-file-base-path": "/home/berg/GitHub/monica-parameters/", # path to monica-parameters
        "path-to-csv-output-dir": "./out/csv-out/"
    },
    "remoteProducer-remoteMonica": {
        "include-file-base-path": "/project/monica-parameters/", # path to monica-parameters
        "path-to-csv-output-dir": "/out/csv-out/"
    }
}


def run():

    config = {
        "mode": "mbm-local-remote",
        "soil_sr": "capnp://login01.cluster.zalf.de:10000",
        "monica_sr": "capnp://login01.cluster.zalf.de:12003",
        "climate_sr": "capnp://login01.cluster.zalf.de:9999",
        "start_row": "0", 
        "end_row": "-1",
        "start_col": "0",
        "end_col": "-1",
        "sim.json": "sim.json",
        "crop.json": "crop.json",
        "site.json": "site.json",
        "setups-file": "sim_setups.csv",
        "run-setups": "[1]",
        "landuse": "data/germany/corine2012_1000_gk5.asc",
        "dem": "data/germany/dem_1000_gk5.asc",
        "slope": "data/germany/slope_1000_gk5.asc",
        "soil": "data/germany/buek200_1000_gk5.asc",
        "bl_gk4": (4450260.1496340110898018,5548602.9990254063159227),
        "tr_gk4": (4578413.4585373029112816,5716481.4597823247313499)
    }
    
    # read commandline args only if script is invoked directly from commandline
    if len(sys.argv) > 1 and __name__ == "__main__":
        for arg in sys.argv[1:]:
            k, v = arg.split("=")
            if k in config:
                config[k] = v

    print("config:", config)

    # select paths 
    paths = PATHS[config["mode"]]

    conMan = common.ConnectionManager()

    soil_service = conMan.connect(config["soil_sr"], soil_data_capnp.Service)
    climate_service = conMan.connect(config["climate_sr"], climate_data_capnp.Service)    
    monica_service = conMan.connect(config["monica_sr"], model_capnp.EnvInstance)    

    # read setup from csv file
    setups = Mrunlib.read_sim_setups(config["setups-file"])
    run_setups = json.loads(config["run-setups"])
    print("read sim setups: ", config["setups-file"])

    #transforms geospatial coordinates from one coordinate reference system to another
    wgs84 = CRS.from_epsg(4326) 
    gk5 = CRS.from_epsg(31469) 
    gk4 = CRS.from_epsg(5678) 
    gk5_to_latlon_transformer = Transformer.from_crs(gk5, wgs84, always_xy=True) 
    gk5_to_gk4 = Transformer.from_crs(gk5, gk4, always_xy=True) 

    # Load grids
    ## note numpy is able to load from a compressed file, ending with .gz or .bz2
    
    # height data for germany
    path_to_dem_grid = config["dem"]
    dem_metadata, _ = Mrunlib.read_header(path_to_dem_grid)
    dem_grid = np.loadtxt(path_to_dem_grid, dtype=int, skiprows=6)
    dem_gk5_interpolate = Mrunlib.create_ascii_grid_interpolator(dem_grid, dem_metadata)
    print("read: ", path_to_dem_grid)
    
    # slope data
    path_to_slope_grid = config["slope"]
    slope_metadata, _ = Mrunlib.read_header(path_to_slope_grid)
    slope_grid = np.loadtxt(path_to_slope_grid, dtype=float, skiprows=6)
    slope_gk5_interpolate = Mrunlib.create_ascii_grid_interpolator(slope_grid, slope_metadata)
    print("read: ", path_to_slope_grid)

    # land use data
    path_to_corine_grid = config["landuse"]
    corine_meta, _ = Mrunlib.read_header(path_to_corine_grid)
    corine_grid = np.loadtxt(path_to_corine_grid, dtype=int, skiprows=6)
    corine_gk5_interpolate = Mrunlib.create_ascii_grid_interpolator(corine_grid, corine_meta)
    print("read: ", path_to_corine_grid)

    # soil data
    path_to_soil_grid = config["soil"]
    soil_metadata, _ = Mrunlib.read_header(path_to_soil_grid)
    soil_grid = np.loadtxt(path_to_soil_grid, dtype=int, skiprows=6)
    print("read: ", path_to_soil_grid)

    sent_env_count = 1
    start_time = time.perf_counter()

    # run calculations for each setup
    for _, setup_id in enumerate(run_setups):

        if setup_id not in setups:
            continue
        start_setup_time = time.perf_counter()      

        setup = setups[setup_id]
        gcm = setup["gcm"]
        rcm = setup["rcm"]
        scenario = setup["scenario"]
        ensmem = setup["ensmem"]
        version = setup["version"]
        crop_id = setup["crop-id"]

        #climate_dataset = climate_service.getDatasetFor({"entries": [{"historical": None}]}).wait().datasets
        climate_dataset = climate_service.getAvailableDatasets().wait().datasets[0].data

        # read template sim.json 
        with open(setup.get("sim.json", config["sim.json"])) as _:
            sim_json = json.load(_)
        # change start and end date acording to setup
        if setup["start_date"]:
            sim_json["climate.csv-options"]["start-date"] = str(setup["start_date"])
        if setup["end_date"]:
            sim_json["climate.csv-options"]["end-date"] = str(setup["end_date"]) 
        sim_json["include-file-base-path"] = paths["include-file-base-path"]

        # read template site.json 
        with open(setup.get("site.json", config["site.json"])) as _:
            site_json = json.load(_)
        # read template crop.json
        with open(setup.get("crop.json", config["crop.json"])) as _:
            crop_json = json.load(_)

        # set the current crop used for this run id
        crop_json["cropRotation"][2] = crop_id

        # create environment template from json templates
        env_template = monica_io3.create_env_json_from_json_config({
            "crop": crop_json,
            "site": site_json,
            "sim": sim_json,
            "climate": ""
        })

        scols = int(soil_metadata["ncols"])
        srows = int(soil_metadata["nrows"])
        scellsize = int(soil_metadata["cellsize"])
        xllcorner = int(soil_metadata["xllcorner"])
        yllcorner = int(soil_metadata["yllcorner"])

        print("All Rows x Cols: " + str(srows) + "x" + str(scols))
        for srow in range(0, srows):
            print(srow,)
            
            if srow < int(config["start_row"]):
                continue
            elif int(config["end_row"]) > 0 and srow > int(config["end_row"]):
                break

            for scol in range(0, scols):
                
                if scol < int(config["start_col"]):
                    continue
                elif int(config["end_col"]) > 0 and scol > int(config["end_col"]):
                    break

                sh_gk5 = yllcorner + (scellsize / 2) + (srows - srow - 1) * scellsize
                sr_gk5 = xllcorner + (scellsize / 2) + scol * scellsize
                
                # check if in bounding box
                sr_gk4, sh_gk4 = gk5_to_gk4.transform(sr_gk5, sh_gk5)
                bl_r, bl_h = config["bl_gk4"]
                tr_r, tr_h = config["tr_gk4"]
                if sr_gk4 < bl_r or sr_gk4 > tr_r or sh_gk4 > tr_h or sh_gk4 < bl_h:
                    continue

                soil_id = int(soil_grid[srow, scol])
                if soil_id == -9999:
                    continue

                height_nn = dem_gk5_interpolate(sr_gk5, sh_gk5)
                slope = slope_gk5_interpolate(sr_gk5, sh_gk5)

                # check if current grid cell is used for agriculture                
                if setup["landcover"]:
                    corine_id = corine_gk5_interpolate(sr_gk5, sh_gk5)
                    if corine_id not in [2,3,4]:
                        continue

                lon, lat = gk5_to_latlon_transformer.transform(sr_gk5, sh_gk5)

                proms = []
                proms.append(soil_service.profilesAt(
                    coord={"lat": lat, "lon": lon},
                    query={"onlyRawData": False, "mandatory": ["soilType", "sand", "clay", "organicCarbon", "bulkDensity"]})
                )
                proms.append(climate_dataset.closestTimeSeriesAt({"lat": lat, "lon": lon}))

                res_list = capnp.join_promises(proms).wait()
                soil_profiles = res_list[0].profiles
                timeseries = res_list[1].timeSeries
                    
                #print("soil:", soil_profile)
                #env_template["params"]["siteParameters"]["SoilProfileParameters"] = soil_profile

                if setup["elevation"]:
                    env_template["params"]["siteParameters"]["heightNN"] = float(height_nn)

                if setup["slope"]:
                    env_template["params"]["siteParameters"]["slope"] = slope / 100.0

                if setup["latitude"]:
                    env_template["params"]["siteParameters"]["Latitude"] = lat

                if setup["CO2"]:
                    env_template["params"]["userEnvironmentParameters"]["AtmosphericCO2"] = float(setup["CO2"])

                if setup["O3"]:
                    env_template["params"]["userEnvironmentParameters"]["AtmosphericO3"] = float(setup["O3"])

                env_template["params"]["simulationParameters"]["UseNMinMineralFertilisingMethod"] = setup["fertilization"]
                env_template["params"]["simulationParameters"]["UseAutomaticIrrigation"] = setup["irrigation"]

                env_template["params"]["simulationParameters"]["NitrogenResponseOn"] = setup["NitrogenResponseOn"]
                env_template["params"]["simulationParameters"]["WaterDeficitResponseOn"] = setup["WaterDeficitResponseOn"]
                env_template["params"]["simulationParameters"]["EmergenceMoistureControlOn"] = setup["EmergenceMoistureControlOn"]
                env_template["params"]["simulationParameters"]["EmergenceFloodingControlOn"] = setup["EmergenceFloodingControlOn"]


                env_template["customId"] = {
                    "setup_id": setup_id,
                    "srow": srow, "scol": scol,
                    "env_id": sent_env_count
                }

                sp = soil_profiles[0]
                env_str = json.dumps(env_template)
                st = common_capnp.StructuredText()
                st.value = env_str
                st.structure.json = None
                monica_res = monica_service.run({"timeSeries": timeseries, "soilProfile": sp, "rest": st }).wait().result

                jres = json.loads(monica_res.as_struct(common_capnp.StructuredText).value)            
                write_csv_file(paths["path-to-csv-output-dir"], jres)

                print("run env ", sent_env_count, " customId: ", env_template["customId"])
                sent_env_count += 1


        stop_setup_time = time.perf_counter()
        print("Setup ", (sent_env_count-1), " envs took ", (stop_setup_time - start_setup_time), " seconds")

    stop_time = time.perf_counter()

    try:
        print("running ", (sent_env_count-1), " envs took ", (stop_time - start_time), " seconds")
        print("exiting run()")
    except Exception:
        raise


def write_csv_file(output_dir, monica_result):

    cid = monica_result.get("customId", {"setup_id": 0, "env_id": 0, "srow": 0, "scol": 0})
    setup_id = str(cid["setup_id"])
    env_id = str(cid["env_id"])
    srow = str(cid["srow"])
    scol = str(cid["scol"])

    path_to_out_dir = output_dir + "/" + setup_id + "/"
    if not os.path.exists(path_to_out_dir):
        try:
            os.makedirs(path_to_out_dir)
        except OSError:
            print("c: Couldn't create dir:", path_to_out_dir, "! Exiting.")
    
    with open(path_to_out_dir + "env-" + env_id + "_srow-" + srow + "_col-" + scol + ".csv", "w", newline='') as _:
        writer = csv.writer(_, delimiter=",")

        for data_ in monica_result.get("data", []):
            results = data_.get("results", [])
            orig_spec = data_.get("origSpec", "")
            output_ids = data_.get("outputIds", [])

            if len(results) > 0:
                writer.writerow([orig_spec.replace("\"", "")])
                for row in monica_io3.write_output_header_rows(output_ids,
                                                            include_header_row=True,
                                                            include_units_row=True,
                                                            include_time_agg=False):
                    writer.writerow(row)

                for row in monica_io3.write_output(output_ids, results):
                    writer.writerow(row)

            writer.writerow([])


def create_output(msg):
    cm_count_to_vals = defaultdict(dict)
    for data in msg.get("data", []):
        results = data.get("results", [])

        is_daily_section = data.get("origSpec", "") == '"daily"'

        for vals in results:
            if "CM-count" in vals:
                cm_count_to_vals[vals["CM-count"]].update(vals)
            elif is_daily_section:
                cm_count_to_vals[vals["Date"]].update(vals)

    cmcs = list(cm_count_to_vals.keys())
    cmcs.sort()
    last_cmc = cmcs[-1]
    if "year" not in cm_count_to_vals[last_cmc]:
        cm_count_to_vals.pop(last_cmc)

    return cm_count_to_vals


def write_row_to_grids(row_col_data, row, ncols, header, path_to_output_dir, path_to_csv_output_dir, setup_id, is_bgr):
    "write grids row by row"
    
    if not hasattr(write_row_to_grids, "nodata_row_count"):
        write_row_to_grids.nodata_row_count = defaultdict(lambda: 0)
        write_row_to_grids.list_of_output_files = defaultdict(list)

    make_dict_nparr = lambda: defaultdict(lambda: np.full((ncols,), -9999, dtype=np.float))
    
    if is_bgr:
        output_grids = {}
        for i in range(1,21):
            output_grids[f'Mois_{i}'] = {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4}
            output_grids[f'STemp_{i}'] = {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4}
        output_keys = ["Mois", "STemp"]
    else:
        output_grids = {
            "sdoy": {"data" : make_dict_nparr(), "cast-to": "int"},
            "ssm03": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "ssm36": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "ssm69": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            
            "s2doy": {"data" : make_dict_nparr(), "cast-to": "int"},
            "s2sm03": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "s2sm36": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "s2sm69": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},

            "sedoy": {"data" : make_dict_nparr(), "cast-to": "int"},
            "sesm03": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "sesm36": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "sesm69": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},

            "s3doy": {"data" : make_dict_nparr(), "cast-to": "int"},
            "s3sm03": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "s3sm36": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "s3sm69": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},

            "s4doy": {"data" : make_dict_nparr(), "cast-to": "int"},
            "s4sm03": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "s4sm36": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "s4sm69": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},

            "s5doy": {"data" : make_dict_nparr(), "cast-to": "int"},
            "s5sm03": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "s5sm36": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "s5sm69": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},

            "s6doy": {"data" : make_dict_nparr(), "cast-to": "int"},
            "s6sm03": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "s6sm36": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "s6sm69": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},

            "s7doy": {"data" : make_dict_nparr(), "cast-to": "int"},
            "s7sm03": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "s7sm36": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "s7sm69": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},

            "hdoy": {"data" : make_dict_nparr(), "cast-to": "int"},
            "hsm03": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "hsm36": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
            "hsm69": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
        }
        output_keys = list(output_grids.keys())

    cmc_to_crop = {}

    is_no_data_row = True
    # skip this part if we write just a nodata line
    if row in row_col_data:
        no_data_cols = 0
        for col in range(0, ncols):
            if col in row_col_data[row]:
                rcd_val = row_col_data[row][col]
                if rcd_val == -9999:
                    no_data_cols += 1
                    continue
                else:
                    cmc_and_year_to_vals = defaultdict(lambda: defaultdict(list))
                    for cell_data in rcd_val:
                        # if we got multiple datasets per cell, iterate over them and aggregate them in the following step
                        for cm_count, data in cell_data.items():
                            for key in output_keys:
                                # store mapping cm_count to crop name for later file name creation
                                if cm_count not in cmc_to_crop and "Crop" in data:
                                    cmc_to_crop[cm_count] = data["Crop"]

                                # only further process/store data we actually received
                                if key in data:
                                    v = data[key]
                                    if isinstance(v, list):
                                        for i, v_ in enumerate(v):
                                            cmc_and_year_to_vals[(cm_count, data["Year"])][f'{key}_{i+1}'].append(v_)
                                    else:
                                        cmc_and_year_to_vals[(cm_count, data["Year"])][key].append(v)
                                # if a key is missing, because that monica event was never raised/reached, create the empty list
                                # so a no-data value is being produced
                                else:
                                    cmc_and_year_to_vals[(cm_count, data["Year"])][key]

                    # potentially aggregate multiple data per cell and finally store them for this row
                    for (cm_count, year), key_to_vals in cmc_and_year_to_vals.items():
                        for key, vals in key_to_vals.items():
                            output_vals = output_grids[key]["data"]
                            if len(vals) > 0:
                                output_vals[(cm_count, year)][col] = sum(vals) / len(vals)
                            else:
                                output_vals[(cm_count, year)][col] = -9999

        is_no_data_row = no_data_cols == ncols

    if is_no_data_row:
        write_row_to_grids.nodata_row_count[setup_id] += 1

    def write_nodata_rows(file_):
        for _ in range(write_row_to_grids.nodata_row_count[setup_id]):
            rowstr = " ".join(["-9999" for __ in range(ncols)])
            file_.write(rowstr +  "\n")

    # iterate over all prepared data for a single row and write row
    for key, y2d_ in output_grids.items():
        y2d = y2d_["data"]
        cast_to = y2d_["cast-to"]
        digits = y2d_.get("digits", 0)
        if cast_to == "int":
            mold = lambda x: str(int(x))
        else:
            mold = lambda x: str(round(x, digits))

        for (cm_count, year), row_arr in y2d.items():
            crop = cmc_to_crop[cm_count] if cm_count in cmc_to_crop else "none"    
            crop = crop.replace("/", "").replace(" ", "")
            path_to_file = path_to_output_dir + crop + "_" + key + "_" + str(year) + "_" + str(cm_count) + ".asc"

            if not os.path.isfile(path_to_file):
                with open(path_to_file, "w") as _:
                    _.write(header)
                    write_row_to_grids.list_of_output_files[setup_id].append(path_to_file)

            with open(path_to_file, "a") as file_:
                write_nodata_rows(file_)
                rowstr = " ".join(["-9999" if int(x) == -9999 else mold(x) for x in row_arr])
                file_.write(rowstr +  "\n")

    # clear the no-data row count when no-data rows have been written before a data row
    if not is_no_data_row:
        write_row_to_grids.nodata_row_count[setup_id] = 0

    # if we're at the end of the output and just empty lines are left, then they won't be written in the
    # above manner because there won't be any rows with data where they could be written before
    # so add no-data rows simply to all files we've written to before
    if is_no_data_row \
        and write_row_to_grids.list_of_output_files[setup_id] \
        and write_row_to_grids.nodata_row_count[setup_id] > 0:
        for path_to_file in write_row_to_grids.list_of_output_files[setup_id]:
            with open(path_to_file, "a") as file_:
                write_nodata_rows(file_)
        write_row_to_grids.nodata_row_count[setup_id] = 0
    
    if row in row_col_data:
        del row_col_data[row]


if __name__ == "__main__":
    run()