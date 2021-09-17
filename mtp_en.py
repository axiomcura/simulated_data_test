# this script calls the R script developed 
# this is a prototype
import subprocess
import numpy as np
import shutil

def parse_output(output_string):
    """ parses the output obtained from MTP-EN script
    
    returns
    -------
    list of tuples containing weights and 
    """

    data = [string_data.split() for string_data in output_string.splitlines()]

    pairs = []
    results = []
    for idx in range(len(data)):
        if len(data[idx]) == 0:
            continue
        elif data[idx][0].startswith("Weight"):
            weight_and_scores = tuple(zip(data[idx], data[idx+1]))
            for weight_score in weight_and_scores:
                pairs.append(weight_score)
        elif data[idx][0].startswith("$"):
            weight_and_score = (data[idx][0], data[idx+1][1:])
            results.append(weight_and_score)
        else:
            continue

    return pairs + results


# def call_mtp_en(n, p1, p2, r1, r2r1ratio, coef1, coef2, beta0, pho_gene, pho_meth, pho_inter, c1, c2):
def call_mtp_en():
    """Calls the MTP-EN algorithm R script 
    
    Arguments
    ---------
    n: int
        Sample size
    p1: int
        total number of features in platform 1
    p2: int
        total number of features in platform 2
    r1: int
        Number of informative features
    r2r1ratio: int
        ratio of informative features from both features 
    coef1: float
        The coeffici9ent of informative features in platform 1
    coef2: float
        The coeffici9ent of informative features in platform 2
    beta0: int
        intercept of the model
    pho_gene: int
        No correalted features
    pho_meth: int
        No correalted features
    pho_interL int
        No correalted features
    c1: int
        Number of informative features that are correlated with platform 1
    c2: int
        Number of informative features that are correlated with platform 2
    Returns
    -------
    """

    # checks if the Rscript executer is installed in the computer
    if not shutil.which("Rscript"):
        raise RuntimeError("Unable to find Rscript executable. Please download R program package")

    # calling MTP-EN script from the dependices folder
    # currently using default paramters of the script (not inputs provded yet)
    cmd_inputs = "Rscript ./dependencies/MTP_EN.r".split() #NOTE: this will change with inputs
    mtp_en = subprocess.run(cmd_inputs, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if mtp_en.returncode != 0:
        print("Error: Error occured within the MTP_EN algorithm")
        print(mtp_en.stderr)
        exit(1)

    mtp_en_output = mtp_en.stdout.decode("utf-8")
    results = parse_output(mtp_en_output)

    #  
    for result in results:
        print(result)

# main script 
if __name__ == "__main__":
    call_mtp_en()