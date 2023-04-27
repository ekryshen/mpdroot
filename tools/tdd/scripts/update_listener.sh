#!/bin/bash

# TDD script
# Written for MPDRoot @JINR Dubna
# First commit Slavomir Hnatic, 4.2023

##########################################################################################################
########################################### Functionality ################################################
###                                                                                                    ###
###                                                                                                    ###
###     DESCRIPTION                                                                                    ###
###                                                                                                    ###
###     The communication between the local machine and the jupyterHub remote server is implemented    ###
###     analogically to method of operation of the signal/slot mechanism (see e.g. Qt library).        ###
###     This script is the "slot" part of the implementation.                                          ###
###     It runs as a listener on the remote server in the background.                                  ###
###                                                                                                    ###
###                                                                                                    ###
###     HOW IT WORKS                                                                                   ###
###                                                                                                    ###
###     -- the script runs only if $LOCAL_LCK_FILE is present                                          ###
###        !!! Important: $LOCAL_LCK_FILE equals $REMOTE_LCK_FILE from eos_rsync.sh !!!                ###
###                                                                                                    ###
###     -- the script executes the update $UPDATE_SCRIPT                                               ###
###                                                                                                    ###
###     -- after cuccessful update the scripts deletes $LOCAL_LCK_FILE                                 ###
###                                                                                                    ###
###                                                                                                    ###
##########################################################################################################

LOCAL_LCK_FILE=update.lck
UPDATE_SCRIPT