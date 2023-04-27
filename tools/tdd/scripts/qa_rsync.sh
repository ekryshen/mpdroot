#!/bin/bash

# TDD script
# Written for MPDRoot @JINR Dubna
# First commit Slavomir Hnatic, 5.2023

##########################################################################################################
########################################### Functionality ################################################
###                                                                                                    ###
###                                                                                                    ###
###     DESCRIPTION                                                                                    ###
###                                                                                                    ###
###     The communication between the local machine and the jupyterHub remote server is implemented    ###
###     analogically to method of operation of the signal/slot mechanism (see e.g. Qt library).        ###
###     This script is the "signal" part of the implementation.                                        ###
###                                                                                                    ###
###     HOW IT WORKS                                                                                   ###
###                                                                                                    ###
###     -- the script checks if (or waits until) the $LCK_FILE is not present                          ###
###        on the remote cluster machine $REMOTE_WORKDIR_ADDRESS to run rsync                          ###
###                                                                                                    ###
###     -- the script rsync's                                                                          ###
###        the list of files $RSYNC_FILES to $REMOTE_WORKDIR_ADDRESS                                   ###
###                                                                                                    ###
###     -- if rsync updated at least 1 file, then 0-byte $REMOTE_LCK_FILE                              ###
###        is created and pushed to the $REMOTE_WORKDIR_ADDRESS                                        ###
###                                                                                                    ###
###                                                                                                    ###
###     NOTE: you should setup automatic login to $REMOTE_SSH_ADDRESS by using ssh key, see e.g.:      ###
###           https://serverfault.com/questions/241588/how-to-automate-ssh-login-with-password         ###
###                                                                                                    ###
###     USE                                                                                            ###
###                                                                                                    ###
###        qa_rsync.sh -u your_login_name                                                              ###
###                                                                                                    ###
###      "your_login_name" is your username on the cluster $REMOTE_SSH_ADDRESS                         ###
###                                                                                                    ###
###                                                                                                    ###
##########################################################################################################

main () {

 # gather login credentials
 while getopts u: flag
 do
    case "${flag}" in
         u) LOGIN_NAME=${OPTARG};;
    esac
 done

 if [ -z "$LOGIN_NAME" ]; then
   echo "Username not set! Exiting..."
   exit 0
 fi 

 # basic variables - files & their addresses
 LCK_FILE=update.lck
 WORKDIR="qa_workdir"
 REMOTE_SSH_ADDRESS="$LOGIN_NAME@hydra.jinr.ru"
 ssh $REMOTE_SSH_ADDRESS "mkdir -p $WORKDIR"
 REMOTE_WORKDIR_ADDRESS="$REMOTE_SSH_ADDRESS:$WORKDIR"
 
 # test for $LCK_FILE presence on the remote machine
 SLEEP_TIME=1 
 SLEEP_TIME_TOTAL=0
 SLEEP_TIME_MAX=60
 test_if_remote_locked
 while [[ "$REMOTE_LOCKED" -eq 1 && "$SLEEP_TIME_TOTAL" -lt "$SLEEP_TIME_MAX" ]]
 do
  echo "Remote machine locked by $LCK_FILE, waiting $SLEEP_TIME sec to unlock..."
  sleep $SLEEP_TIME
  ((++SLEEP_TIME_TOTAL))
  test_if_remote_locked
 done

 if [[ "$SLEEP_TIME_TOTAL" -ge "$SLEEP_TIME_MAX" ]]; then
  echo "Remote machine not updating! Exiting..."
  exit 11
 fi

 # remote server not write locked, running rsync
 run_rsync

}

function test_if_remote_locked () {
 ssh $REMOTE_SSH_ADDRESS "test -e $WORKDIR/$LCK_FILE" && REMOTE_LOCKED=1 || REMOTE_LOCKED=0
}

function run_rsync () {
 RSYNC_COMMAND=$(rsync -ai -e ssh *.root $REMOTE_WORKDIR_ADDRESS)
 if [ -n "${RSYNC_COMMAND}" ]; then
  echo "Changes were made by rsync, locking the remote machine for update & exiting."
  ssh $REMOTE_SSH_ADDRESS "touch $WORKDIR/$LCK_FILE"
 else
   echo "No changes were made by rsync on the remote machine."
 fi
}

main "$@"
