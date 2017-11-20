import re

def Capture_Cromwell_ID():
    File = open("/projects/mgc/Project_1/ram/CromwellWDL_WorkFlow_Development/WorkflowCodes/Genomics_MGC_GenomeGPS_CromwelWDL/qsubs/output.txt", "r")

    rows = [line.strip() for line in File]

    # test = rows[9]

    for i in range(len(rows)):
        searchObj = re.search(r'"id.*"(.*)"', rows[i], re.M | re.I)
        if searchObj:
            Cromwell_ID = searchObj.group(1)
    File.close()
    return Cromwell_ID



def main():
    x = Capture_Cromwell_ID()
    print x
main()
