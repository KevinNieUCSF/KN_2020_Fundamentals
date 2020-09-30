
#import modules
import sys
import numpy as np
import scipy
import math
import pypdb as pd
import collections
import matplotlib.pyplot as plt
import pandas as pds
from pylab import rcParams

"""this code is made to recreate figure 1a and 1b from Laitaoja, et al. "Zinc Coordination Spheres in Protein Structures" """

def parsepdb(iqmain, iqart): #general pipeline, tentative
    """fetch IDs from pdb that contain zinc, returns two list of IDs, one for zinc and one for zinc x-ray artifacts"""
    qmain=pd.make_query(iqmain, querytype='AdvancedKeywordQuery')
    qart=pd.make_query(iqart, querytype= 'AdvancedKeywordQuery')
    qmainl=[]
    qartfinal=[]
    for hit in pd.do_search(qmain):
        qmainl.append(hit)
    for hit in pd.do_search(qart):
        qartfinal.append(hit)
    return qmainl,qartfinal
def pdblistfilter(qmainl,iqnmr,iqxray):
    """filters out zinc hits for NMR or x-ray methodology"""
    qnmr=pd.make_query(iqnmr, querytype='ExpTypeQuery')
    qxray=pd.make_query(iqxray, querytype='ExpTypeQuery')
    qnmrl=[]
    qxrayl=[]
    for hit in pd.do_search(qnmr):
        qnmrl.append(hit)
    for hit in pd.do_search(qxray):
        qxrayl.append(hit)
    qnmrfinal=[]
    qxrayfinal=[]
    for element in qmainl:
        if element in qnmrl:
            qnmrfinal.append(element)
        if element in qxrayl:
            qxrayfinal.append(element)
    return qnmrfinal, qxrayfinal
def molcounternmr(test,nmr):
    try:
        if (float(test['polymer']['@weight']))<=2500:
            nmr[2.5]+=1
        elif (float(test['polymer']['@weight']))<=5000:
            nmr[5]+=1
        elif (float(test['polymer']['@weight']))<=7500:
            nmr[7.5]+=1
        elif (float(test['polymer']['@weight']))<=10000:
            nmr[10]+=1
        elif (float(test['polymer']['@weight']))<=12500:
            nmr[12.5]+=1
        elif (float(test['polymer']['@weight']))<=15000:
            nmr[15]+=1
        elif (float(test['polymer']['@weight']))<=17500:
            nmr[17.5]+=1
        elif (float(test['polymer']['@weight']))<=20000:
            nmr[20]+=1
        elif (float(test['polymer']['@weight']))<=22500:
            nmr[22.5]+=1
        elif (float(test['polymer']['@weight']))<=25000:
            nmr[25]+=1
        elif (float(test['polymer']['@weight']))<=27500:
            nmr[27.5]+=1
        elif (float(test['polymer']['@weight']))<=30000:
            nmr[30]+=1
        elif (float(test['polymer']['@weight']))<=32500:
            nmr[32.5]+=1
        elif (float(test['polymer']['@weight']))<=35000:
            nmr[35]+=1
        elif (float(test['polymer']['@weight']))<=37500:
            nmr[37.5]+=1
    except:
        pass
    try:
        polymer=test['polymer']
        mol=0
        if type(polymer) is list:
            for poly in polymer:
                mol+=float(poly['@weight'])
        if mol<=2500:
            nmr[2.5]+=1
        elif mol<=5000:
            nmr[5]+=1
        elif mol<=7500:
            nmr[7.5]+=1
        elif mol<=10000:
            nmr[10]+=1
        elif mol<=12500:
            nmr[12.5]+=1
        elif mol<=15000:
            nmr[15]+=1
        elif mol<=17500:
            nmr[17.5]+=1
        elif mol<=20000:
            nmr[20]+=1
        elif mol<=22500:
            nmr[22.5]+=1
        elif mol<=25000:
            nmr[25]+=1
        elif mol<=27500:
            nmr[27.5]+=1
        elif mol<=30000:
            nmr[30]+=1
        elif mol<=32500:
            nmr[32.5]+=1
        elif mol<=35000:
            nmr[35]+=1
        elif mol<=37500:
            nmr[37.5]+=1
    except:
        pass
    else:
        pass
    return nmr
def molcounterxray(test,xray):
    try:
        if (float(test['polymer']['@weight']))<=20000:
            xray[20]+=1
        elif (float(test['polymer']['@weight']))<=40000:
            xray[40]+=1
        elif (float(test['polymer']['@weight']))<=60000:
            xray[60]+=1
        elif (float(test['polymer']['@weight']))<=80000:
            xray[80]+=1
        elif (float(test['polymer']['@weight']))<=100000:
            xray[100]+=1
        elif (float(test['polymer']['@weight']))<=120000:
            xray[120]+=1
        elif (float(test['polymer']['@weight']))<=140000:
            xray[140]+=1
        elif (float(test['polymer']['@weight']))<=160000:
            xray[160]+=1
        elif (float(test['polymer']['@weight']))<=180000:
            xray[180]+=1
        elif (float(test['polymer']['@weight']))<=200000:
            xray[200]+=1
        elif (float(test['polymer']['@weight']))<=220000:
            xray[220]+=1
        elif (float(test['polymer']['@weight']))<=240000:
            xray[240]+=1
        elif (float(test['polymer']['@weight']))<=260000:
            xray[260]+=1
        elif (float(test['polymer']['@weight']))<=280000:
            xray[280]+=1
        elif (float(test['polymer']['@weight']))<=1000000:
            xray[1000]+=1
    except:
        pass
    try:
        polymer=test['polymer']
        mol=0
        if type(polymer) is list:
            for poly in polymer:
                mol+=float(poly['@weight'])
        if mol<=20000:
            xray[20]+=1
        elif mol<=40000:
            xray[40]+=1
        elif mol<=60000:
            xray[60]+=1
        elif mol<=80000:
            xray[80]+=1
        elif mol<=100000:
            xray[100]+=1
        elif mol<=120000:
            xray[120]+=1
        elif mol<=140000:
            xray[140]+=1
        elif mol<=160000:
            xray[160]+=1
        elif mol<=180000:
            xray[180]+=1
        elif mol<=200000:
            xray[200]+=1
        elif mol<=220000:
            xray[220]+=1
        elif mol<=240000:
            xray[240]+=1
        elif mol<=260000:
            xray[260]+=1
        elif mol<=280000:
            xray[280]+=1
        elif mol<=1000000:
            xray[1000]+=1
    except:
        pass
    else:
        pass
    return xray
def molcounterart(test,art):
    try:
        if (float(test['polymer']['@weight']))<=20000:
            art[20]+=1
        elif (float(test['polymer']['@weight']))<=40000:
            art[40]+=1
        elif (float(test['polymer']['@weight']))<=60000:
            art[60]+=1
        elif (float(test['polymer']['@weight']))<=80000:
            art[80]+=1
        elif (float(test['polymer']['@weight']))<=100000:
            art[100]+=1
        elif (float(test['polymer']['@weight']))<=120000:
            art[120]+=1
        elif (float(test['polymer']['@weight']))<=140000:
            art[140]+=1
        elif (float(test['polymer']['@weight']))<=160000:
            art[160]+=1
        elif (float(test['polymer']['@weight']))<=180000:
            art[180]+=1
        elif (float(test['polymer']['@weight']))<=200000:
            art[200]+=1
        elif (float(test['polymer']['@weight']))<=220000:
            art[220]+=1
        elif (float(test['polymer']['@weight']))<=240000:
            art[240]+=1
        elif (float(test['polymer']['@weight']))<=260000:
            art[260]+=1
        elif (float(test['polymer']['@weight']))<=280000:
            art[280]+=1
        elif (float(test['polymer']['@weight']))<=1000000:
            art[1000]+=1
    except:
        pass
    try:
        polymer=test['polymer']
        mol=0
        if type(polymer) is list:
            for poly in polymer:
                mol+=float(poly['@weight'])
        if mol<=20000:
            art[20]+=1
        elif mol<=40000:
            art[40]+=1
        elif mol<=60000:
            art[60]+=1
        elif mol<=80000:
            art[80]+=1
        elif mol<=100000:
            art[100]+=1
        elif mol<=120000:
            art[120]+=1
        elif mol<=140000:
            art[140]+=1
        elif mol<=160000:
            art[160]+=1
        elif mol<=180000:
            art[180]+=1
        elif mol<=200000:
            art[200]+=1
        elif mol<=220000:
            art[220]+=1
        elif mol<=240000:
            art[240]+=1
        elif mol<=260000:
            art[260]+=1
        elif mol<=280000:
            art[280]+=1
        elif mol<=1000000:
            art[1000]+=1
    except:
        pass
    else:
        pass
    return art
def fetchmol(qnmrfinal,qxrayfinal,qartfinal):
    """gets the molecular weight of each hit for each generated ID list and then generates a count
    for the coming bar graph. mol count info stored in dicts. The info of what weight each ID is is not recorded,
    rather, their existence is recorded if they are within a range of kDa"""
    nmrmolrange=np.linspace(2.5,37.5,15,endpoint=True)
    nmr={i:0 for i in nmrmolrange}
    xrayartmolrange=np.linspace(20,280,14,endpoint=True)
    xray={i:0 for i in xrayartmolrange}
    xray[1000]=0
    art={i:0 for i in xrayartmolrange}
    art[1000]=0
    for hit in qnmrfinal:
        try:
            test=pd.get_all_info(hit)
            nmr=molcounternmr(test,nmr)
        except:
            continue
    for hit in qxrayfinal:
        try:
            test=pd.get_all_info(hit)
            xray=molcounterxray(test,xray)
        except:
            continue
    for hit in qartfinal:
        try:
            test=pd.get_all_info(hit)
            art=molcounterart(test,art)
        except:
            continue
    print("NMR Mol Weight Tally:"+str(nmr))
    print("Xray Mol Weight Tally:"+str(xray))
    print("Artifact Mol Weight Tally:"+str(art))
    return nmr, xray, art
def graph(nmr,xray,art):
    """generate a graph using the dictionaries generated from parsepdb"""
    #nmr_dataframe = pds.dataframe.from_dict(nmr)
    #xray_dataframe = pds.dataframe.from_dict(xray)
    #artifact_dataframe = pds.dataframe.from_dict(art)

    #test data:
    NMR ={2.5: 1021, 5.0: 265, 7.5: 207, 10.0: 245, 12.5: 139, 15.0: 92, 17.5: 70, 20.0: 59, 22.5: 18, 25.0: 13, 27.5: 7, 30.0: 14, 32.5: 4, 35.0: 0, 37.5: 0}
    xray = {20.0: 15227, 40.0: 6710, 60.0: 4030, 80.0: 1022, 100.0: 498, 120.0: 395, 140.0: 79, 160.0: 43, 180.0: 49, 200.0: 47, 220.0: 60, 240.0: 11, 260.0: 3, 280.0: 6, 1000: 434}
    art = {20.0: 459, 40.0: 192, 60.0: 121, 80.0: 22, 100.0: 21, 120.0: 22, 140.0: 0, 160.0: 0, 180.0: 1, 200.0: 8, 220.0: 0, 240.0: 0, 260.0: 0, 280.0: 0, 1000: 0}


    plt.rcParams.update({'font.size': 22})

    nmr_plot = plt.bar(list(NMR.keys()), NMR.values(), color='b', tick_label=list(NMR.keys()))
    plt.rcParams["figure.figsize"] = (30,10)
    plt.xlabel("Structure Molecular Mass (kDa)")
    plt.ylabel("Number of NMR Structures")
    plt.title('NMR')
    nmr_plot.show()

    plt.bar(list(xray.keys()), xray.values(), color='b', width = 5, tick_label=list(xray.keys()))
    shifted = [x + 5 for x in list(art.keys())]
    plt.bar(shifted, art.values(), color='r', width = 5)
    plt.rcParams["figure.figsize"] = (22,10)
    plt.xlabel("Structure Molecular Mass (kDa)")
    plt.ylabel("Number of Xray Structures")
    plt.title('Xray')
        #plt.xlim(0, 290)
    #plt.ylim(0, 7500)
    plt.show()

    return None
def main():
    iqmain="zinc"
    iqnmr="NMR"
    iqxray="X-RAY"
    iqart="zinc X-RAY artifact"
    print('Querying PDB...')
    qmainl,qartfinal=parsepdb(iqmain, iqart)
    print('Filtering Lists...')
    qnmrfinal,qxrayfinal=pdblistfilter(qmainl,iqnmr,iqxray)
    print("Fetching molweights...")
    nmr,xray,art=fetchmol(qnmrfinal,qxrayfinal,qartfinal)
    graph(nmr,xray,art)

main()

"""nmrmolrange=np.linspace(2.5,37.5,15,endpoint=True)
nmr={i:0 for i in nmrmolrange}
xrayartmolrange=np.linspace(20,280,14,endpoint=True)
xray={i:0 for i in xrayartmolrange}
xray[1000]=0
art={i:0 for i in xrayartmolrange}
art[1000]=0
hit='1pvn'
test=pd.get_all_info(hit)
print(test)
HAHAHAHAA THIS DOESNT WORK FOR SOME PDB IDS HOOHOOHEEEHEE I LOVE CODING"""
#print(test['molDescription']['structureId']['polymer']['@weight'])



"""if dict polymer is list
def get_weight(polymer)
    return polymer['@weight']

def total_weight(pdb):
polymer = pdb[polymer]
total_weight = 0.0
if type(polymer) is list:
    for  poly in polymer:
        total_weight = total_weight + get_weight(poly)
else:
    total_weight = get_weight(polymer)


 """
