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

import time
import os
import math
import json
import csv
import copy
#from io import StringIO
from datetime import date, timedelta
from collections import defaultdict
import sys
import zmq

import sqlite3
import sqlite3 as cas_sq3
import numpy as np
from pyproj import CRS, transform

import monica_io3
import soil_io3
import monica_run_lib as Mrunlib

PATHS = {
    # adjust the local path to your environment
    "mbm-local-remote": {
        "include-file-base-path": "/home/berg/GitHub/monica-parameters/", # path to monica-parameters
        "path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/dwd_core_ensemble/", # mounted path to archive or hard drive with climate data 
        "monica-path-to-climate-dir": "/monica_data/climate-data/dwd_core_ensemble/csvs_dwd_core_ensemble/", # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./monica-data/data/", # mounted path to archive or hard drive with data 
        "path-to-projects-dir": "./monica-data/data/projects/",
        "path-debug-write-folder": "./debug-out/",
    },
    "remoteProducer-remoteMonica": {
        "include-file-base-path": "/project/monica-parameters/", # path to monica-parameters
        "path-to-climate-dir": "/data/dwd_core_ensemble/", # mounted path to archive or hard drive with climate data 
        "monica-path-to-climate-dir": "/monica_data/climate-data/dwd_core_ensemble/csvs_dwd_core_ensemble/", # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./monica-data/data/", # mounted path to archive or hard drive with data 
        "path-to-projects-dir": "./monica-data/data/projects/", # mounted path to archive or hard drive with project data 
        "path-debug-write-folder": "/out/debug-out/",
    }
}

DEFAULT_HOST = "login01.cluster.zalf.de" # "localhost" #
DEFAULT_PORT = "6666"
RUN_SETUP = "[3]"
SETUP_FILE = "sim_setups.csv"
PROJECT_FOLDER = "monica-germany/"
DATA_SOIL_DB = "germany/buek200.sqlite"
DATA_GRID_HEIGHT = "germany/dem_1000_gk5.asc" 
DATA_GRID_SLOPE = "germany/slope_1000_gk5.asc"
DATA_GRID_LAND_USE = "germany/landuse_1000_gk5.asc"
DATA_GRID_SOIL = "germany/BUEK200_1000_gk5.asc"
TEMPLATE_PATH_LATLON = "{path_to_climate_dir}/latlon-to-rowcol.json"
TEMPLATE_PATH_CLIMATE_CSV = "{gcm}/{rcm}/{scenario}/{ensmem}/{version}/row-{crow}/col-{ccol}.csv"
GEO_TARGET_GRID=31469 #proj4 -> 3-degree gauss-kruger zone 5 (=Germany) https://epsg.io/5835 ###https://epsg.io/31469

DEBUG_DONOT_SEND = False
DEBUG_WRITE = False
DEBUG_ROWS = 10
DEBUG_WRITE_FOLDER = "./debug_out"
DEBUG_WRITE_CLIMATE = False

# some values in these templates will be overwritten by the setup 
TEMPLATE_SIM_JSON="sim.json" 
TEMPLATE_CROP_JSON="crop.json"
TEMPLATE_SITE_JSON="site.json"

# commandline parameters e.g "server=localhost port=6666 shared_id=2"
def run_producer(server = {"server": None, "port": None}, shared_id = None):
    "main"

    context = zmq.Context()
    socket = context.socket(zmq.PUSH) # pylint: disable=no-member
    #config_and_no_data_socket = context.socket(zmq.PUSH)

    config = {
        "mode": "mbm-local-remote",
        "server-port": server["port"] if server["port"] else DEFAULT_PORT,
        "server": server["server"] if server["server"] else DEFAULT_HOST,
        "start-row": "0", 
        "end-row": "-1",
        "sim.json": TEMPLATE_SIM_JSON,
        "crop.json": TEMPLATE_CROP_JSON,
        "site.json": TEMPLATE_SITE_JSON,
        "setups-file": SETUP_FILE,
        "run-setups": RUN_SETUP,
        "shared_id": shared_id
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
    # open soil db connection
    soil_db_con = sqlite3.connect(paths["path-to-data-dir"] + DATA_SOIL_DB)
    #soil_db_con = cas_sq3.connect(paths["path-to-data-dir"] + DATA_SOIL_DB) #CAS.
    # connect to monica proxy (if local, it will try to connect to a locally started monica)
    socket.connect("tcp://" + config["server"] + ":" + str(config["server-port"]))

    # read setup from csv file
    setups = Mrunlib.read_sim_setups(config["setups-file"])
    run_setups = json.loads(config["run-setups"])
    print("read sim setups: ", config["setups-file"])

    #transforms geospatial coordinates from one coordinate reference system to another
    # transform wgs84 into gk5
    wgs84 = CRS.from_epsg(4326) #proj4 -> (World Geodetic System 1984 https://epsg.io/4326)
    gk5 = CRS.from_epsg(GEO_TARGET_GRID) 

    # Load grids
    ## note numpy is able to load from a compressed file, ending with .gz or .bz2
    
    # height data for germany
    path_to_dem_grid = paths["path-to-data-dir"] + DATA_GRID_HEIGHT 
    dem_metadata, _ = Mrunlib.read_header(path_to_dem_grid)
    dem_grid = np.loadtxt(path_to_dem_grid, dtype=int, skiprows=6)
    dem_gk5_interpolate = Mrunlib.create_ascii_grid_interpolator(dem_grid, dem_metadata)
    print("read: ", path_to_dem_grid)
    
    # slope data
    path_to_slope_grid = paths["path-to-data-dir"] + DATA_GRID_SLOPE
    slope_metadata, _ = Mrunlib.read_header(path_to_slope_grid)
    slope_grid = np.loadtxt(path_to_slope_grid, dtype=float, skiprows=6)
    slope_gk5_interpolate = Mrunlib.create_ascii_grid_interpolator(slope_grid, slope_metadata)
    print("read: ", path_to_slope_grid)

    # land use data
    path_to_corine_grid = paths["path-to-data-dir"] + DATA_GRID_LAND_USE
    corine_meta, _ = Mrunlib.read_header(path_to_corine_grid)
    corine_grid = np.loadtxt(path_to_corine_grid, dtype=int, skiprows=6)
    corine_gk5_interpolate = Mrunlib.create_ascii_grid_interpolator(corine_grid, corine_meta)
    print("read: ", path_to_corine_grid)

    # soil data
    path_to_soil_grid = paths["path-to-data-dir"] + DATA_GRID_SOIL
    soil_metadata, _ = Mrunlib.read_header(path_to_soil_grid)
    soil_grid = np.loadtxt(path_to_soil_grid, dtype=int, skiprows=6)
    print("read: ", path_to_soil_grid)

    cdict = {}
    # path to latlon-to-rowcol.json
    path = TEMPLATE_PATH_LATLON.format(path_to_climate_dir=paths["path-to-climate-dir"])
    climate_data_gk5_interpolator = Mrunlib.create_climate_geoGrid_interpolator_from_json_file(path, wgs84, gk5, cdict)
    print("created climate_data to gk5 interpolator: ", path)

    sent_env_count = 1
    start_time = time.perf_counter()

    listOfClimateFiles = set()
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

        # read template sim.json 
        with open(setup.get("sim.json", config["sim.json"])) as _:
            sim_json = json.load(_)
        # change start and end date acording to setup
        if setup["start_date"]:
            sim_json["climate.csv-options"]["start-date"] = str(setup["start_date"])
        if setup["end_date"]:
            sim_json["climate.csv-options"]["end-date"] = str(setup["end_date"]) 
        sim_json["include-file-base-path"] = paths["include-file-base-path"]

        if setup["bgr"]:
            if setup["nc_mode"]:
                sim_json["output"]["events"] = sim_json["output"]["nc-bgr-events"]
            else:
                sim_json["output"]["events"] = sim_json["output"]["bgr-events"]
        else:
            if setup["nc_mode"]:
                sim_json["output"]["events"] = sim_json["output"]["nc-events"]

        sim_json["output"]["obj-outputs?"] = not setup["nc_mode"] and not setup["bgr"]

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

        # set shared id in template
        if config["shared_id"]:
            env_template["sharedId"] = config["shared_id"]

        scols = int(soil_metadata["ncols"])
        srows = int(soil_metadata["nrows"])
        scellsize = int(soil_metadata["cellsize"])
        xllcorner = int(soil_metadata["xllcorner"])
        yllcorner = int(soil_metadata["yllcorner"])

        #unknown_soil_ids = set()
        soil_id_cache = {}
        print("All Rows x Cols: " + str(srows) + "x" + str(scols))
        for srow in range(0, srows):
            print(srow,)
            
            if srow < int(config["start-row"]):
                continue
            elif int(config["end-row"]) > 0 and srow > int(config["end-row"]):
                break

            for scol in range(0, scols):
                soil_id = int(soil_grid[srow, scol])
                if soil_id == -9999:
                    continue

                if soil_id in soil_id_cache:
                    soil_profile = soil_id_cache[soil_id]
                else:
                    soil_profile = soil_io3.soil_parameters(soil_db_con, soil_id)
                    soil_id_cache[soil_id] = soil_profile

                if len(soil_profile) == 0:
                    print("row/col:", srow, "/", scol, "has unknown soil_id:", soil_id)
                    #unknown_soil_ids.add(soil_id)
                    continue
                
                #get coordinate of clostest climate element of real soil-cell
                sh_gk5 = yllcorner + (scellsize / 2) + (srows - srow - 1) * scellsize
                sr_gk5 = xllcorner + (scellsize / 2) + scol * scellsize
                #inter = crow/ccol encoded into integer
                crow, ccol = climate_data_gk5_interpolator(sr_gk5, sh_gk5)

                # check if current grid cell is used for agriculture                
                if setup["landcover"]:
                    corine_id = corine_gk5_interpolate(sr_gk5, sh_gk5)
                    if corine_id not in [2,3,4]:
                        continue

                height_nn = dem_gk5_interpolate(sr_gk5, sh_gk5)
                slope = slope_gk5_interpolate(sr_gk5, sh_gk5)

                env_template["params"]["userCropParameters"]["__enable_T_response_leaf_expansion__"] = setup["LeafExtensionModifier"]

                    
                #print("soil:", soil_profile)
                env_template["params"]["siteParameters"]["SoilProfileParameters"] = soil_profile

                # setting groundwater level
                if setup["groundwater-level"]:
                    groundwaterlevel = 20
                    layer_depth = 0
                    for layer in soil_profile:
                        if layer.get("is_in_groundwater", False):
                            groundwaterlevel = layer_depth
                            #print("setting groundwaterlevel of soil_id:", str(soil_id), "to", groundwaterlevel, "m")
                            break
                        layer_depth += Mrunlib.get_value(layer["Thickness"])
                    env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepthMonth"] = 3
                    env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepth"] = [max(0, groundwaterlevel - 0.2) , "m"]
                    env_template["params"]["userEnvironmentParameters"]["MaxGroundwaterDepth"] = [groundwaterlevel + 0.2, "m"]
                    
                # setting impenetrable layer
                if setup["impenetrable-layer"]:
                    impenetrable_layer_depth = Mrunlib.get_value(env_template["params"]["userEnvironmentParameters"]["LeachingDepth"])
                    layer_depth = 0
                    for layer in soil_profile:
                        if layer.get("is_impenetrable", False):
                            impenetrable_layer_depth = layer_depth
                            #print("setting leaching depth of soil_id:", str(soil_id), "to", impenetrable_layer_depth, "m")
                            break
                        layer_depth += Mrunlib.get_value(layer["Thickness"])
                    env_template["params"]["userEnvironmentParameters"]["LeachingDepth"] = [impenetrable_layer_depth, "m"]
                    env_template["params"]["siteParameters"]["ImpenetrableLayerDepth"] = [impenetrable_layer_depth, "m"]

                if setup["elevation"]:
                    env_template["params"]["siteParameters"]["heightNN"] = float(height_nn)

                if setup["slope"]:
                    env_template["params"]["siteParameters"]["slope"] = slope / 100.0

                if setup["latitude"]:
                    clat, _ = cdict[(crow, ccol)]
                    env_template["params"]["siteParameters"]["Latitude"] = clat

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

                env_template["csvViaHeaderOptions"] = sim_json["climate.csv-options"]
                
                subpath_to_csv = TEMPLATE_PATH_CLIMATE_CSV.format(gcm=gcm, rcm=rcm, scenario=scenario, ensmem=ensmem, version=version, crow=str(crow), ccol=str(ccol))
                env_template["pathToClimateCSV"] = [paths["monica-path-to-climate-dir"] + subpath_to_csv]
                if setup["incl_hist"]:
                    hist_subpath_to_csv = TEMPLATE_PATH_CLIMATE_CSV.format(gcm=gcm, rcm=rcm, scenario="historical", ensmem=ensmem, version=version, crow=str(crow), ccol=str(ccol))
                    env_template["pathToClimateCSV"].insert(0, paths["monica-path-to-climate-dir"] + hist_subpath_to_csv)
                print(env_template["pathToClimateCSV"])
                if DEBUG_WRITE_CLIMATE :
                    listOfClimateFiles.add(subpath_to_csv)

                env_template["customId"] = {
                    "setup_id": setup_id,
                    "srow": srow, "scol": scol,
                    "crow": int(crow), "ccol": int(ccol),
                    "soil_id": soil_id,
                    "bgr": setup["bgr"],
                    "env_id": sent_env_count
                }

                if not DEBUG_DONOT_SEND :
                    socket.send_json(env_template)
                    print("sent env ", sent_env_count, " customId: ", env_template["customId"])

                sent_env_count += 1

                # write debug output, as json file
                if DEBUG_WRITE:
                    debug_write_folder = paths["path-debug-write-folder"]
                    if not os.path.exists(debug_write_folder):
                        os.makedirs(debug_write_folder)
                    if sent_env_count < DEBUG_ROWS  :

                        path_to_debug_file = debug_write_folder + "/row_" + str(sent_env_count-1) + "_" + str(setup_id) + ".json" 

                        if not os.path.isfile(path_to_debug_file):
                            with open(path_to_debug_file, "w") as _ :
                                _.write(json.dumps(env_template))
                        else:
                            print("WARNING: Row ", (sent_env_count-1), " already exists")
            #print("unknown_soil_ids:", unknown_soil_ids)

            #print("crows/cols:", crows_cols)
        stop_setup_time = time.perf_counter()
        print("Setup ", (sent_env_count-1), " envs took ", (stop_setup_time - start_setup_time), " seconds")

    stop_time = time.perf_counter()

    # write summary of used json files
    if DEBUG_WRITE_CLIMATE:
        debug_write_folder = paths["path-debug-write-folder"]
        if not os.path.exists(debug_write_folder):
            os.makedirs(debug_write_folder)

        path_to_climate_summary = debug_write_folder + "/climate_file_list" + ".csv"
        with open(path_to_climate_summary, "w") as _:
            _.write('\n'.join(listOfClimateFiles))

    try:
        print("sending ", (sent_env_count-1), " envs took ", (stop_time - start_time), " seconds")
        #print("ran from ", start, "/", row_cols[start], " to ", end, "/", row_cols[end]
        print("exiting run_producer()")
    except Exception:
        raise

if __name__ == "__main__":
    run_producer()