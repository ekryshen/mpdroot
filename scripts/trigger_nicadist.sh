#! /bin/bash
# Script to install mpdroot based on tag
# Written for MPDRoot @JINR Dubna
# First commit Jan Busa, 12.2021
if [[ $CI_COMMIT_REF_NAME == dev ]]; then
  curl --request POST --header "PRIVATE-TOKEN: $NICADIST_CI_TOKEN" \
  --header "Content-Type: application/json" \
  --data '{ "ref": "master", "variables": [ {"key": "NICADIST_PACKAGE", "value": "mpdroot-dev-tmp"} ] }' \
  "https://git.jinr.ru/api/v4/projects/982/pipeline"
elif [[ $CI_COMMIT_REF_NAME == *-tmp ]]; then
  curl --request POST --header "PRIVATE-TOKEN: $NICADIST_CI_TOKEN" \
  --header "Content-Type: application/json" \
  --data '{ "ref": "master", "variables": [ {"key": "NICADIST_PACKAGE", "value": "mpdroot-'$CI_COMMIT_REF_NAME'"} ] }' \
  "https://git.jinr.ru/api/v4/projects/982/pipeline"
else
  NICADIST_TAG=mpdroot-$CI_COMMIT_REF_NAME
  curl --request POST --header "PRIVATE-TOKEN: $NICADIST_CI_TOKEN" \
  --header "Content-Type: application/json" \
  --data '{ "ref": "master", "tag_name": "'$NICADIST_TAG'" }' \
  "https://git.jinr.ru/api/v4/projects/982/repository/tags"
fi;
