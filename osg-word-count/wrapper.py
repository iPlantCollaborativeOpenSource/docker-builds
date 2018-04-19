#!/usr/bin/env python

from __future__ import print_function
from os.path import expanduser, basename

import csv
import itertools
import json
import os
import requests
import socket
import subprocess
import sys

# A clean way to print to stderr.
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# A class to load configuration settings.
class Config:
    def __init__(self, config_json):
        self.irods_host = self.extract_setting(config_json, "irods_host")
        self.irods_port = self.extract_setting(config_json, "irods_port", default=1247)
        self.irods_job_user = self.extract_setting(config_json, "irods_job_user")
        self.input_ticket_list = self.extract_setting(config_json, "input_ticket_list", default="input_ticket.list")
        self.output_ticket_list = self.extract_setting(config_json, "output_ticket_list", default="output_ticket.list")
        self.status_update_url = self.extract_setting(config_json, "status_update_url")
        self.stdout = self.extract_setting(config_json, "stdout", default="out.txt")
        self.stderr = self.extract_setting(config_json, "stderr", default="err.txt")

    def extract_setting(self, config_json, attribute, default=None):
        if attribute in config_json:
            return config_json[attribute]
        elif default is not None:
            return default
        else:
            eprint("required configuration setting {0} not provided".format(attribute))
            sys.exit(1)

# This is a simple context manager class designed to make it easier to read and write iRODS ticket list files.
class TicketListReader:
    def __init__(self, path):
        self.path = path

    def __enter__(self):
        self.fd = open(self.path, "rb")
        self.r = csv.reader(itertools.ifilter(lambda l: l[0] != '#', self.fd))
        return self.r

    def __exit__(self, type, value, traceback):
        self.fd.close()

# This is a simple class used to send job status update messages.
class JobStatusUpdater:
    def __init__(self, url):
        self.url = url

    def print_update(self, status, message):
        print("{0}: {1}".format(status, message))

    def send_update(self, status, message):
        self.print_update(status, message)
        body = {"status": status, "message": message, "hostname": socket.gethostname()}
        r = requests.post(self.url, json=body)
        if r.status_code < 200 or r.status_code > 299:
            eprint("unable to send job status update: {0} {1}".format(status, message))

    def failed(self, message):
        self.send_update("failed", message)

    def running(self, message):
        self.send_update("running", message)

    def completed(self, message):
        self.send_update("completed", message)

# Load the configuration file.
def load_config():
    with open("config.json", "r") as f:
        return Config(json.load(f))

# Ensure that a directory exists.
def ensuredir(path, mode):
    if not os.path.exists(path):
        os.mkdir(path, mode)

# Initialize iRODS.
def init_irods(host, port):
    irods_config = {
        "irods_user_name": "anonymous",
        "irods_host": host,
        "irods_port": port,
        "irods_zone_name": ""
    }
    ensuredir(expanduser("~/.irods"), 0755)
    with open(expanduser("~/.irods/irods_environment.json"), "w") as f:
        json.dump(irods_config, f, indent=4, separators=(",", ": "))

# Download a file or directory from iRODS.
def download_file(ticket, src):
    rc = subprocess.call(["iget", "-rt", ticket, src])
    if rc != 0:
        raise Exception("could not download {0}".format(src))

# Download a set of files referenced in a ticket list file from iRODS, returning a list of downloaded files.
def download_files(ticket_list_path):
    result = []
    with TicketListReader(ticket_list_path) as tlr:
        for ticket, src in tlr:
            download_file(ticket, src)
            result.append(basename(src))
    return result

# Run the word count job.
def run_job(input_files, output_filename, error_filename):
    with open(output_filename, "w") as out, open(error_filename, "w") as err:
        rc = subprocess.call(["wc"] + input_files, stdout=out, stderr=err)
        if rc != 0:
            raise Exception("wc returned exit code {0}".format(rc))

# Upload a file or directory to iRODS.
def upload_file(ticket, owner, src, dest):
    rc = subprocess.call(["iput", "-rt", ticket, src, dest])
    if rc != 0:
        raise Exception("could not upload {0} to {1}".format(src, dest))
    fullpath = os.path.join(dest, src)
    rc = subprocess.call(["ichmod", "own", owner, fullpath])
    if rc != 0:
        raise Exception("could not change the owner of {0} to {1}".format(fullpath, owner))
    rc = subprocess.call(["ichmod", "null", "anonymous", fullpath])
    if rc != 0:
        raise Exception("could not remove anonymous permissions on {0}".format(fullpath))

# Upload a set of files to the directories referenced in a ticket list file to iRODS.
def upload_files(ticket_list_path, owner, files):
    with TicketListReader(ticket_list_path) as tlr:
        for ticket, dest in tlr:
            for src in files:
                upload_file(ticket, owner, src, dest)

if __name__ == "__main__":
    # Load the configuration file.
    cfg = load_config()

    # Create a job status updater and send the running update.
    updater = JobStatusUpdater(cfg.status_update_url)
    updater.running("configuration successfully loaded")

    # Initialize iRODS.
    try:
        updater.running("initializing the iRODS connection")
        init_irods(cfg.irods_host, cfg.irods_port)
    except Exception as e:
        updater.failed("unable to initialize the iRODS connection: {0}".format(e))
        sys.exit(1)

    # Download the files from iRODS, assuming that all tickets refer to regular files for now.
    try:
        updater.running("downloading the input files")
        input_files = download_files(cfg.input_ticket_list)
    except Exception as e:
        updater.failed("unable to download input files: {0}".format(e))
        sys.exit(1)

    # Keep track of whether or not the job has failed.
    job_failed = False

    # Run the job.
    try:
        updater.running("processing the input files")
        run_job(input_files, cfg.stdout, cfg.stderr)
    except Exception as e:
        updater.running("job encountered an error: {0}".format(e))
        job_failed = True

    # Upload the output file.
    try:
        updater.running("uploading the output files")
        upload_files(cfg.output_ticket_list, cfg.irods_job_user, [cfg.stdout, cfg.stderr])
    except Exception as e:
        updater.running("unable to upload output file: {0}".format(e))
        job_failed = True

    # Send the final job status update.
    if job_failed:
        updater.failed("job failed; see prior status update messages for details")
        sys.exit(1)
    else:
        updater.completed("job completed successfully")