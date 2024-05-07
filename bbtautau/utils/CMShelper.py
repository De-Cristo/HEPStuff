import os, sys
import subprocess
import ROOT
import argparse
import time

def dasgoclient(query):
    cmd = 'dasgoclient -query="file dataset={0}"'.format(query)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print("Reading {0} files".format(len(result.stdout.strip().split('\n'))))
    print(result.stdout.strip().split('\n'))
    return result.stdout.strip().split('\n')

def remote_reading(database):
    redirector = ["root://cms-xrd-global.cern.ch/", "root://xrootd-cms.infn.it/", "root://cmsxrootd.fnal.gov/"]
    
    # the exeption method is not really working for some reason
    try:
        infile_list = [redirector[0] + data for data in database]
        df = ROOT.RDataFrame('Events', infile_list)
    except Exception as e:
        print("An error occurred while creating the RDataFrame object:\n", e)
        print(redirector[0]+" reading failed, try: "+redirector[1])
        try:
            infile_list = [redirector[1] + data for data in database]
            df = ROOT.RDataFrame('Events', infile_list)
        except Exception as e:
            print("An error occurred while creating the RDataFrame object:\n", e)
            print(redirector[1]+" reading failed, try: "+redirector[2])
            try:
                infile_list = [redirector[2] + data for data in database]
                df = ROOT.RDataFrame('Events', infile_list)
            except Exception as e:
                print("An error occurred while creating the RDataFrame object:\n", e)
                sys.exit(1)
                
    return df

