#! /bin/bash

nohup python auto_golden_gate_server.py 2>&1 1> server_log.txt &
disown
 
