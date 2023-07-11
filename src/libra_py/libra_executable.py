#***********************************************************
# * Copyright (C) 2021 Brendan Smith and Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

import sys
import json
import logging

import libra_py.recipes as recipes


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler('libra_executable.log')
file_handler.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


def read_libra_jobs_json(input_file):
    """ This function reads the libra_jobs.json file and returns a python dictionary

    Args:
        input_file (string): The name of the .json file that contains the user-defined
            parameters for using the nbra/workflows in a black-box manner

    Returns:
        params ( dictionary ): Parameters for carrying out the nbra/workflows in a 
            black-box manner
    """

    logger.debug("Entered into the read_libra_jobs_json function of libra_executables.py")

    with open(input_file) as f:
        file_lines = f.readlines()
        f.close()
 
    json_string = ""
    for file_line in file_lines:
        file_line = file_line.replace('\n', '')
        json_string = json_string + file_line

    params = json.loads(json_string)

    logger.debug("Debugging the json -> python params dict")
    logger.debug(params)

    return params 



def main(json_parameters_filename):
    """ Call this function to run your libra 'executable'

    Args:
        json_parameters_filename (string): The name of the .json file 
            that contains the user-defined parameters for executing a 
            Libra capability

    Returns:
        None: but carrys out the Libra capability in a 
            black-box manner
    """

    logger.debug("Entered into the main() function of libra_executables.py")

    params = read_libra_jobs_json('libra_jobs.json')

    jobs = params["jobs"]
    for job in jobs:

        logger.debug(F"job {job}\n")
        if job["job_type"] == "step2":
            logger.debug("Doing some step2 job")

            params_indx = job["job_params"]

            logger.debug("Generating step2 submit_template.slm")
            recipes.step2_recipes.make_step2_submit_template(params["step2_params"][params_indx])

            logger.debug("Initializing step2 job folders and files")
            recipes.step2_recipes.run_step2_jobs(params["step2_params"][params_indx])

        if job["job_type"] == "pdos":
            logger.debug("Doing some pdos job")

            params_indx = job["job_params"]

            logger.debug(params["pdos_params"][params_indx])
            recipes.pdos_recipes.compute_pdos(params["pdos_params"][params_indx])


if __name__ == '__main__':
    main("libra_jobs.json")
