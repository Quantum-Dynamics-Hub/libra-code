import sys
import json
import logging

from libra_py.workflows.nbra import recipes


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

    with open(input_file) as f:
        file_lines = f.readlines()
        f.close()
 
    json_string = ""
    for file_line in file_lines:
        file_line = file_line.replace('\n', '')
        json_string = json_string + file_line

    params = json.loads(json_string)
    logger.debug(params)

    return params 




#================ Do the job =============
params = read_libra_jobs_json('libra_jobs.json')

jobs = params["jobs"]
for job in jobs:

    logger.debug(F"job {job}\n")
    if job["job_type"] == "step2":
        logger.debug("Doing step2")
        params_indx = job["job_params"]
        logger.debug("Generating step2 submit_template.slm")
        recipes.generate_step2_submit_template(params["step2_params"][params_indx])
        logger.debug("Initializing step2 job folders and files")
        recipes.initialize_step2_jobs(params["step2_params"][params_indx])

    #if job["job_type"] == "pdos":
        #logger.debug("Doing pDOS\n")
        #params_indx = job["job_params"]
        #logger.debug(params["pdos_params"][params_indx])
        #recipes.pdos.run(job["pdos_params"][params_indx])   
