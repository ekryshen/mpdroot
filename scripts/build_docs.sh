#!/bin/bash

# Build gitlab pages for mpdroot project
# Written for MPDRoot @JINR Dubna
# First commit Jan Busa Jr., 04.2022

######################################################################
######################### Functionality ##############################
###                                                                ###
###  Scripts build mpdroot based documentation as well as MkDocs   ###
###  It works in 3 steps:                                          ###
###   1) using sources from tools/documentation/doxygen            ###
###      a html version of reference manual is created             ###
###   2) using sources from tools/documentation/webpage            ###
###      a mkdocs version of reference manual is created           ###
###   3) combining files from 1) and 2) directory pages is created ###
###      and published in git.jinr.ru pages                        ###
###  This documentation can be then found at page                  ###
###             http://nica.pages.jinr.ru/mpdroot/                 ###
###                                                                ###
######################################################################

# build mkdocs documentation and place result into folder pages
rm -rf public/
cd tools/documentation/webpage
mkdocs build
mv site/ ../../../public
cd ../../..
cp public/index.html public/404.html
## doxygen build is disabled for now
# copy doxygen configuration files
##cp tools/documentation/doxygen/buildDoxy.sh .
##cp tools/documentation/doxygen/MPDrootDoxy .
### build doxygen documentation
##./buildDoxy.sh
##mv html/ public/refman
##rm -f public/refman/*.md5 # clean unnecessary files
##rm -f public/refman/*.map # clean unnecessary files
##rm -f buildDoxy.sh
##rm -f MPDrootDoxy
