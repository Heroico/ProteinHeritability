__author__ = 'heroico'

import re

KEY_REML_VG = "V(G)"
KEY_REML_VG_SE = "SE(V(G))"
KEY_REML_VE = "V(e)"
KEY_REML_VE_SE = "SE(V(e))"
KEY_REML_VP = "Vp"
KEY_REML_VP_SE = "SE(Vp)"
KEY_REML_VG_TO_VP = "V(G)_TO_Vp"
KEY_REML_VG_TO_VP_SE = "SE(V(G)_TO_Vp)"
KEY_REML_LOGL = "logL"
KEY_REML_LOGL0 = "logL0"
KEY_REML_LRT = "LRT"
KEY_REML_DF = "df"
KEY_REML_PVAL = "Pval"
KEY_REML_N = "n"

vg_re = re.compile(r"[-+]?\d*\.\d+|\d+")
def LoadREMLFromFile(file_name):
    result = {}
    with open(file_name, 'rb') as file:
        lines = []
        for line in file:
            lines.append(line)

        vg_results = vg_re.findall(lines[1])
        result[KEY_REML_VG] = vg_results[0]
        result[KEY_REML_VG_SE] = vg_results[1]

        ve_results = vg_re.findall(lines[2])
        result[KEY_REML_VE] = ve_results[0]
        result[KEY_REML_VE_SE] = ve_results[1]

        vp_results = vg_re.findall(lines[3])
        result[KEY_REML_VP] = vp_results[0]
        result[KEY_REML_VP_SE] = vp_results[1]

        vg_to_vp_results = vg_re.findall(lines[4])
        result[KEY_REML_VG_TO_VP] = vg_to_vp_results[0]
        result[KEY_REML_VG_TO_VP_SE] = vg_to_vp_results[1]

        lg_results = vg_re.findall(lines[5])
        result[KEY_REML_LOGL] = lg_results[0]

        lg0_results = vg_re.findall(lines[6])
        result[KEY_REML_LOGL0] = lg0_results[1]

        lrt_results = vg_re.findall(lines[7])
        result[KEY_REML_LRT] = lrt_results[0]

        df_results = vg_re.findall(lines[8])
        result[KEY_REML_DF] = df_results[0]

        pval_results = vg_re.findall(lines[9])
        result[KEY_REML_PVAL] = pval_results[0]

        n_results = vg_re.findall(lines[10])
        result[KEY_REML_N] = n_results[0]
    return result

def NameFromFileName(fileName):
    name = fileName[20:-4]
    return name

def OrderFromName(name):
    results = vg_re.findall(name)
    result = results[0]
    order = int(result)
    return order

KEY_RESULT_NAME = "name"
KEY_RESULT_REML = "reml"
KEY_RESULT_ORDER = "order"
def BuildResultsFromFiles(files):
    results = []
    for file_name in files:
        reml = LoadREMLFromFile(file_name)
        name = NameFromFileName(file_name)
        order = OrderFromName(name)
        result = {KEY_RESULT_NAME: name, KEY_RESULT_REML:reml, KEY_RESULT_ORDER:order}
        results.append(result)
    results=sorted(results, key=lambda result: result[KEY_RESULT_ORDER])
    return results

RESULTS_HEADER=header = "PHENO,V(G),SE(V(G)),V(e),SE(V(e),Vp,SE(vp),V(G)_to_Vp,SE(V(G)_to_Vp),logL,logL0,LRT,df,Pval,n\n"
def PrintResultsToFile(results,file_name='Out/reml_results.csv'):
    with open(file_name, 'w+') as out:
        out.write(RESULTS_HEADER)
        for result in results:
            data=result[KEY_RESULT_REML]
            name=result[KEY_RESULT_NAME]
            fields = []
            fields.append(name)
            fields.append(data[KEY_REML_VG])
            fields.append(data[KEY_REML_VG_SE])
            fields.append(data[KEY_REML_VE])
            fields.append(data[KEY_REML_VE_SE])
            fields.append(data[KEY_REML_VP])
            fields.append(data[KEY_REML_VP_SE])
            fields.append(data[KEY_REML_VG_TO_VP])
            fields.append(data[KEY_REML_VG_TO_VP_SE])
            fields.append(data[KEY_REML_LOGL])
            fields.append(data[KEY_REML_LOGL0])
            fields.append(data[KEY_REML_LRT])
            fields.append(data[KEY_REML_DF])
            fields.append(data[KEY_REML_PVAL])
            fields.append(data[KEY_REML_N])
            line = ",".join(fields)+"\n"
            out.write(line)


if __name__ == "__main__":
    import glob
    import argparse

    parser = argparse.ArgumentParser(description='Figure out REML from phenos.')
    parser.add_argument("--file_input_prefix",
                        help="pattern of reml files")
    parser.add_argument("--reml_output",
                        help="name of output file")
    args = parser.parse_args()

    pattern = './Intermediate/reml*'
    if args.file_input_prefix and len(args.file_input_prefix):
        pattern = args.file_input_prefix
    files = glob.glob(pattern)
    results = BuildResultsFromFiles(files)

    fine_name = 'Out/reml_results.csv'
    if args.reml_output and len(args.reml_output):
        file_name = args.reml_output
    PrintResultsToFile(results, file_name)