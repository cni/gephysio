#!/usr/bin/env python

import os
import json

# Parse a config file
def parse_config(config_json_file):
    """
    Take a json file, read and return the contents
    """

    if not os.path.isfile(config_json_file):
        raise ValueError('No config file could be found!')

    # Read the config json file
    with open(config_json_file, 'r') as jsonfile:
        config = json.load(jsonfile)

    return config

if __name__ == '__main__':

    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--config_file', type=str, dest="config_file", default='/flywheel/v0/config.json', help='Full path to the input json config file.')
    ap.add_argument('--debug', '-d', type=bool, dest="debug", default=False, help='show debug print statements')
    args = ap.parse_args()

    config = parse_config(args.config_file)
    flywheel_input = '/flywheel/v0/input/'
    flywheel_output = '/flywheel/v0/output'

    import flywheel

    #import flywheel.GearContext as context
    #context.init_logging()
    #context.log_config()

    fw = flywheel.Client(config['inputs']['api-key']['key'])
    acqID = config['inputs']['physio']['hierarchy']['id']
    acquisition = fw.get_acquisition(acqID)

    for f in acquisition['files']:
        if f['type'] == 'dicom':
            nvols = f.info['NumberOfTemporalPositions']
            TR = f.info['RepetitionTime']     # in milliseconds
            break

    # Establish time base for the PPG data
    #niftiinfo = config['inputs']['nifti']['object']['info']
    #TR = niftiinfo['RepititionTime']
    #nvols = niftiinfo['NumberOfTemporalPositions']

    # Flag for generating physio data plot
    plot_flag = 1 if config['config']['Snapshot'] else 0
    plot_start = config['config']['start_time']
    plot_window = config['config']['window_length']

    # Unzip physio data 
    physiofile = config['inputs']['physio']['location']['path']
    print("physiofile: %s\n" % physiofile)
    os.system("unzip %s -d %s" % (physiofile, flywheel_input))
    physio_dir = os.path.join(flywheel_input, os.path.splitext(os.path.basename(physiofile))[0])
    print("physio_dir: %s\n" % physio_dir)
    
    cmd = ("run_gephysio.sh /opt/mcr/v95 %s %s %d %d %d %d %d" %( physio_dir, flywheel_output, TR, nvols, plot_flag, plot_start, plot_window))
    print(cmd)
    os.system(cmd)


