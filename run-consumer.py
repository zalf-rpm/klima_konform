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

from collections import defaultdict, OrderedDict
import csv
from datetime import datetime
import json
import gc
import numpy as np
import os
from pyproj import CRS, Transformer
import shutil
import sqlite3
import sys
import timeit
import types
import zmq

import monica_io3
import soil_io3
import monica_run_lib as Mrunlib

PATHS = {
    "mbm-local-remote": {
        "path-to-data-dir": "data/",
        "path-to-output-dir": "out/",
        "path-to-csv-output-dir": "csv-out/"
    },
    "remoteConsumer-remoteMonica": {
        "path-to-data-dir": "data/",
        "path-to-final-output-dir": "/out/out/",
        "path-to-output-dir": "/scratch/klima_konform/out/out/",
        #"path-to-csv-output-dir": "/out/csv-out/",
        "path-to-csv-output-dir": "/scratch/klima_konform/out/csv-out/"
    }
}
TEMPLATE_SOIL_PATH = "{local_path_to_data_dir}germany/buek200_1000_25832_etrs89-utm32n.asc"
TEMPLATE_LANDUSE_PATH = "{local_path_to_data_dir}germany/landuse_1000_31469_gk5.asc"
USE_LANDUSE = False

def create_output(msg):
    indexval_and_section_to_vals = defaultdict(dict)
    for data in msg.get("data", []):
        results = data.get("results", [])

        is_crop_section = data.get("origSpec", "") == '"crop"'
        is_daily_section = data.get("origSpec", "") == '"daily"'
        is_yearly_section = data.get("origSpec", "") == '"yearly"'

        for vals in results:
            if is_crop_section and "CM-count" in vals and "Harvest-year" in vals:
                indexval_and_section_to_vals[(f'{vals["CM-count"]}_Harvest-year_{vals["Harvest-year"]}', "season")].update(vals)
            elif is_daily_section and "Date" in vals:
                indexval_and_section_to_vals[(vals["Date"], "daily")].update(vals)
            elif is_yearly_section and "Year" in vals:
                indexval_and_section_to_vals[(vals["Year"], "yearly")].update(vals)
            elif "Index" in vals:
                indexval_and_section_to_vals[(vals["Index"], "undef")].update(vals)

    #cmcs = list([t[0] for t in cm_count_to_vals.keys()])
    #cmcs.sort()
    #last_cmc = cmcs[-1]
    #if "Year" not in cm_count_to_vals[last_cmc]:
    #    cm_count_to_vals.pop(last_cmc)

    return indexval_and_section_to_vals


def write_row_to_grids(row_col_data, row, ncols, header, path_to_output_dir, path_to_csv_output_dir, setup_id):
    "write grids row by row"
    
    if not hasattr(write_row_to_grids, "nodata_row_count"):
        write_row_to_grids.nodata_row_count = defaultdict(lambda: 0)
        write_row_to_grids.list_of_output_files = defaultdict(list)

    make_dict_nparr = lambda: defaultdict(lambda: np.full((ncols,), -9999, dtype=np.float))
    
    output_grids = {
        ("Yield", "season"): {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        ("Groundwater_recharge_sum_season", "season"): {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        ("Groundwater_recharge_sum_year", "yearly"): {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        ("AET_sum_season", "season"): {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        ("AET_sum_year", "yearly"): {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        ("SOC_0-30_avg_avg_year", "yearly"): {"data" : make_dict_nparr(), "cast-to": "float", "digits": 4},
        ("Season_length", "season"): {"data" : make_dict_nparr(), "cast-to": "int"},
        ("LAI_max_season", "season"): {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        #"TraDef": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 2},
        #"HeatRed": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 2},
        #"FrostRed": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 2},
    }
    output_keys = list(output_grids.keys())

    indexval_and_section_to_crop = {}

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
                    indexval_and_section_to_vals = defaultdict(lambda: defaultdict(list))
                    for cell_data in rcd_val:
                        # if we got multiple datasets per cell, iterate over them and aggregate them in the following step
                        for (indexval, section), data in cell_data.items():
                            for key in [key for (key, sec) in output_keys if sec == section]:
                                # store mapping cm_count to crop name for later file name creation
                                if (indexval, section) not in indexval_and_section_to_crop and "Crop" in data:
                                    indexval_and_section_to_crop[(indexval, section)] = data["Crop"]

                                # only further process/store data we actually received
                                if key in data:
                                    v = data[key]
                                    if isinstance(v, list):
                                        for i, v_ in enumerate(v):
                                            indexval_and_section_to_vals[(indexval, section)][f'{key}_{i+1}'].append(v_)
                                    else:
                                        indexval_and_section_to_vals[(indexval, section)][key].append(v)
                                # if a key is missing, because that monica event was never raised/reached, create the empty list
                                # so a no-data value is being produced
                                else:
                                    indexval_and_section_to_vals[(indexval, section)][key]

                    # potentially aggregate multiple data per cell and finally store them for this row
                    for (indexval, section), key_to_vals in indexval_and_section_to_vals.items():
                        for key, vals in key_to_vals.items():
                            output_vals = output_grids[(key, section)]["data"]
                            if len(vals) > 0:
                                output_vals[indexval][col] = sum(vals) / len(vals)
                            else:
                                output_vals[indexval][col] = -9999

        is_no_data_row = no_data_cols == ncols

    if is_no_data_row:
        write_row_to_grids.nodata_row_count[setup_id] += 1

    def write_nodata_rows(file_):
        for _ in range(write_row_to_grids.nodata_row_count[setup_id]):
            rowstr = " ".join(["-9999" for __ in range(ncols)])
            file_.write(rowstr +  "\n")

    # iterate over all prepared data for a single row and write row
    for (key, section), y2d_ in output_grids.items():
        y2d = y2d_["data"]
        cast_to = y2d_["cast-to"]
        digits = y2d_.get("digits", 0)
        if cast_to == "int":
            mold = lambda x: str(int(x))
        else:
            mold = lambda x: str(round(x, digits))

        for indexval, row_arr in y2d.items():
            crop = indexval_and_section_to_crop[(indexval, section)] + "_" if (indexval, section) in indexval_and_section_to_crop else ""    
            crop = crop.replace("/", "").replace(" ", "")
            path_to_file = path_to_output_dir + crop + key + "_" + str(indexval) + ".asc"

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


def run_consumer(leave_after_finished_run = True, server = {"server": None, "port": None}, shared_id = None):
    "collect data from workers"

    config = {
        "mode": "mbm-local-remote",
        "port": server["port"] if server["port"] else "7777",
        "server": server["server"] if server["server"] else "login01.cluster.zalf.de", 
        "gk4-tl-br-bounds": "[[4450260.1496340110898018,5716481.4597823247313499],[4578413.4585373029112816,5548602.9990254063159227]]",
        "shared_id": shared_id,
        "timeout": 600000 # 10 minutes
    }

    if len(sys.argv) > 1 and __name__ == "__main__":
        for arg in sys.argv[1:]:
            k,v = arg.split("=")
            if k in config:
                config[k] = v

    paths = PATHS[config["mode"]]

    if not "out" in config:
        config["out"] = paths["path-to-output-dir"]
    if not "csv-out" in config:
        config["csv-out"] = paths["path-to-csv-output-dir"]

    print("consumer config:", config)

    context = zmq.Context()
    if config["shared_id"]:
        socket = context.socket(zmq.DEALER)
        socket.setsockopt(zmq.IDENTITY, config["shared_id"])
    else:
        socket = context.socket(zmq.PULL)

    socket.connect("tcp://" + config["server"] + ":" + config["port"])
    socket.RCVTIMEO = config["timeout"]
    leave = False
    write_normal_output_files = True

    path_to_soil_grid = TEMPLATE_SOIL_PATH.format(local_path_to_data_dir=paths["path-to-data-dir"])
    soil_epsg_code = int(path_to_soil_grid.split("/")[-1].split("_")[2])
    soil_crs = CRS.from_epsg(soil_epsg_code)
    soil_metadata, header = Mrunlib.read_header(path_to_soil_grid)
    soil_grid_template = np.loadtxt(path_to_soil_grid, dtype=int, skiprows=6)

    scols = int(soil_metadata["ncols"])
    srows = int(soil_metadata["nrows"])
    scellsize = int(soil_metadata["cellsize"])
    xllcorner = int(soil_metadata["xllcorner"])
    yllcorner = int(soil_metadata["yllcorner"])
    snodata_value = int(soil_metadata["nodata_value"])

    if config["gk4-tl-br-bounds"] is None or config["gk4-tl-br-bounds"] == "":
        bounds = None
    else:
        gk4 = CRS.from_epsg(5678)
        soil_crs_to_gk4_transformer = Transformer.from_crs(soil_crs, gk4, always_xy=True) 
        ((tl_r_gk4, tl_h_gk4), (br_r_gk4, br_h_gk4)) = json.loads(config["gk4-tl-br-bounds"])

        bounds = {"tl": [srows, scols], "br": [0, 0]}
        #bounds_gk4 = {"tl": [srows, scols], "br": [0, 0]}

        for srow in range(0, srows):
            #print(srow,)
            
            sh_upper_col0_soil_crs = yllcorner + scellsize + (srows - srow - 1) * scellsize
            sh_lower_col0_soil_crs = yllcorner + (srows - srow - 1) * scellsize
            sr_col0_soil_crs = xllcorner + (scellsize / 2)
            row_upper_r_gk4, row_upper_h_gk4 = soil_crs_to_gk4_transformer.transform(sr_col0_soil_crs, sh_upper_col0_soil_crs)
            row_lower_r_gk4, row_lower_h_gk4 = soil_crs_to_gk4_transformer.transform(sr_col0_soil_crs, sh_lower_col0_soil_crs)

            if row_lower_h_gk4 < br_h_gk4 or row_upper_h_gk4 > tl_h_gk4:
                continue 

            for scol in range(0, scols):

                sh_col_soil_crs = yllcorner + (scellsize / 2) + (srows - srow - 1) * scellsize
                sr_left_col_soil_crs = xllcorner + scol * scellsize
                sr_right_col_soil_crs = xllcorner + scellsize + scol * scellsize
                col_left_r_gk4, col_left_h_gk4 = soil_crs_to_gk4_transformer.transform(sr_left_col_soil_crs, sh_col_soil_crs)
                col_right_r_gk4, col_right_h_gk4 = soil_crs_to_gk4_transformer.transform(sr_right_col_soil_crs, sh_col_soil_crs)

                if col_left_r_gk4 < tl_r_gk4 or col_right_r_gk4 > br_r_gk4:
                    #print(col_left_r_gk4, "<", tl_r_gk4, "or", col_right_r_gk4, ">", br_r_gk4)
                    continue 

                if bounds["tl"][0] > srow:
                    bounds["tl"][0] = srow
                    #bounds_gk4["tl"][0] = row_upper_h_gk4
                if bounds["tl"][1] > scol:
                    bounds["tl"][1] = scol
                    #bounds_gk4["tl"][1] = col_left_r_gk4
                if bounds["br"][0] < srow:
                    bounds["br"][0] = srow
                    #bounds_gk4["br"][0] = row_lower_h_gk4
                if bounds["br"][1] < scol:
                    bounds["br"][1] = scol
                    #bounds_gk4["br"][1] = col_right_r_gk4

        print("bounds:", bounds)
        #print(bounds_gk4)
        soil_grid_template= soil_grid_template[bounds["tl"][0]:bounds["br"][0]+1, bounds["tl"][1]:bounds["br"][1]+1]
    
        header = "ncols " + str(soil_grid_template.shape[1]) + "\n" + \
            "nrows " + str(soil_grid_template.shape[0]) + "\n" + \
            "xllcorner " + str(int(soil_metadata["xllcorner"]) + scellsize * bounds["tl"][1]) + "\n" + \
            "yllcorner " + str(int(soil_metadata["yllcorner"]) + scellsize * (srows - bounds["br"][0])) + "\n" + \
            "cellsize " + str(scellsize) + "\n" + \
            "nodata_value " + str(snodata_value) + "\n" 


    if USE_LANDUSE:
        path_to_landuse_grid = TEMPLATE_LANDUSE_PATH.format(local_path_to_data_dir=paths["path-to-data-dir"])
        landuse_epsg_code = int(path_to_landuse_grid.split("/")[-1].split("_")[2])
        landuse_crs = CRS.from_epsg(landuse_epsg_code)
        landuse_transformer = Transformer.from_crs(soil_crs, landuse_crs)
        landuse_meta, _ = Mrunlib.read_header(path_to_landuse_grid)
        landuse_grid = np.loadtxt(path_to_landuse_grid, dtype=int, skiprows=6)
        landuse_interpolate = Mrunlib.create_ascii_grid_interpolator(landuse_grid, landuse_meta)

        for srow in range(0, srows):
            #print(srow)
            for scol in range(0, scols):
                soil_id = soil_grid_template[srow, scol]
                if soil_id == -9999:
                    continue

                #get coordinate of clostest climate element of real soil-cell
                sh = yllcorner + (scellsize / 2) + (srows - srow - 1) * scellsize
                sr = xllcorner + (scellsize / 2) + scol * scellsize

                # check if current grid cell is used for agriculture                
                lur, luh = landuse_transformer(sh, sr)
                landuse_id = landuse_interpolate(lur, luh)
                if landuse_id not in [2,3,4]:
                    soil_grid_template[srow, scol] = -9999

        print("filtered through LANDUSE")

    #set all data values to one, to count them later
    soil_grid_template[soil_grid_template != -9999] = 1
    #set all no-data values to 0, to ignore them while counting
    soil_grid_template[soil_grid_template == -9999] = 0

    #count cols in rows
    datacells_per_row = np.sum(soil_grid_template, axis=1)

    setup_id_to_data = defaultdict(lambda: {
        "end_row": -1,
        "nrows": soil_grid_template.shape[0],
        "ncols": soil_grid_template.shape[1],
        "header": header,
        "out_dir_exists": False,
        "row-col-data": defaultdict(lambda: defaultdict(list)),
        "datacell-count": datacells_per_row.copy(),
        "next-row": 0
    })

    def process_message(msg):
        if len(msg["errors"]) > 0:
            print("There were errors in message:", msg, "\nSkipping message!")
            return


        if not hasattr(process_message, "wnof_count"):
            process_message.wnof_count = 0
            process_message.setup_count = 0

        leave = False

        if not write_normal_output_files:
            custom_id = msg["customId"]
            setup_id = custom_id["setup_id"]
            is_nodata = custom_id["nodata"]

            data = setup_id_to_data[setup_id]

            row = custom_id["row"] 
            col = custom_id["col"]

            debug_msg = "received work result " + str(process_message.received_env_count) + " customId: " + str(msg.get("customId", "")) \
            + " next row: " + str(data["next-row"]) \
            + " cols@row to go: " + str(data["datacell-count"][row]) + "@" + str(row) + " cells_per_row: " + str(datacells_per_row[row])#\
            #+ " rows unwritten: " + str(data["row-col-data"].keys()) 
            print(debug_msg)
            #debug_file.write(debug_msg + "\n")
            if not is_nodata:
                data["row-col-data"][row][col].append(create_output(msg))
            data["datacell-count"][row] -= 1

            process_message.received_env_count = process_message.received_env_count + 1

            #while data["next-row"] in data["row-col-data"] and data["datacell-count"][data["next-row"]] == 0:
            while data["datacell-count"][data["next-row"]] == 0:
                
                path_to_out_dir = config["out"] + str(setup_id) + "/"
                path_to_csv_out_dir = config["csv-out"] + str(setup_id) + "/"
                if not data["out_dir_exists"]:
                    if os.path.isdir(path_to_out_dir) and os.path.exists(path_to_out_dir):
                        data["out_dir_exists"] = True
                    else:
                        try:
                            os.makedirs(path_to_out_dir)
                            data["out_dir_exists"] = True
                        except OSError:
                            print("c: Couldn't create dir:", path_to_out_dir, "! Exiting.")
                            exit(1)
                    if os.path.isdir(path_to_csv_out_dir) and os.path.exists(path_to_csv_out_dir):
                        data["out_dir_exists"] = True
                    else:
                        try:
                            os.makedirs(path_to_csv_out_dir)
                            data["out_dir_exists"] = True
                        except OSError:
                            print("c: Couldn't create dir:", path_to_csv_out_dir, "! Exiting.")
                            exit(1)
                
                write_row_to_grids(data["row-col-data"], data["next-row"], data["ncols"], data["header"], \
                    path_to_out_dir, path_to_csv_out_dir, setup_id)
                
                debug_msg = "wrote row: "  + str(data["next-row"]) + " next-row: " + str(data["next-row"]+1) + " rows unwritten: " + str(list(data["row-col-data"].keys()))
                print(debug_msg)
                #debug_file.write(debug_msg + "\n")
                
                data["next-row"] += 1 # move to next row (to be written)

                if leave_after_finished_run \
                and (data["next-row"] > data["nrows"]-1 or data["next-row"] > data["end_row"]): 
                    process_message.setup_count += 1
                
        elif write_normal_output_files:

            if msg.get("type", "") in ["jobs-per-cell", "no-data", "setup_data"]:
                #print "ignoring", result.get("type", "")
                return

            print("received work result ", process_message.received_env_count, " customId: ", str(msg.get("customId", "")))

            custom_id = msg["customId"]
            setup_id = custom_id["setup_id"]
            row = custom_id["row"]
            col = custom_id["col"]
            
            process_message.wnof_count += 1

            path_to_out_dir = config["out"] + str(setup_id) + "/" + str(row) + "/"
            if not os.path.exists(path_to_out_dir):
                try:
                    os.makedirs(path_to_out_dir)
                except OSError:
                    print("c: Couldn't create dir:", path_to_out_dir, "! Exiting.")
                    exit(1)

            #with open("out/out-" + str(i) + ".csv", 'wb') as _:
            with open(path_to_out_dir + "col-" + str(col) + ".csv", "w", newline='') as _:
                writer = csv.writer(_, delimiter=",")

                for data_ in msg.get("data", []):
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

            process_message.received_env_count = process_message.received_env_count + 1

        return leave

    process_message.received_env_count = 1

    while not leave:
        try:
            #start_time_recv = timeit.default_timer()
            msg = socket.recv_json(encoding="latin-1")
            #elapsed = timeit.default_timer() - start_time_recv
            #print("time to receive message" + str(elapsed))
            #start_time_proc = timeit.default_timer()
            leave = process_message(msg)
            #elapsed = timeit.default_timer() - start_time_proc
            #print("time to process message" + str(elapsed))
        except zmq.error.Again as _e:
            print('no response from the server (with "timeout"=%d ms) ' % socket.RCVTIMEO)

            # copy directory to output directory
            if write_normal_output_files:
                if not os.path.exists(paths["path-to-final-output-dir"]):
                    try:
                        os.makedirs(paths["path-to-final-output-dir"])
                    except OSError:
                        print("c: Couldn't create dir:", paths["path-to-final-output-dir"], "! Exiting.")
                        exit(1)
                shutil.copytree(config["out"], paths["path-to-final-output-dir"])
            return
        except Exception as e:
            print("Exception:", e)
            #continue

    print("exiting run_consumer()")
    #debug_file.close()

if __name__ == "__main__":
    run_consumer()


