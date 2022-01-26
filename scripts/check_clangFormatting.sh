#!/bin/bash

# Clang formatting check for Gitlab CI pipeline
# Written for MPDRoot @JINR Dubna
# First commit Slavomir Hnatic, 11.2021

######################################################################
######################### Functionality ##############################
###                                                                ###
###  1.If any code from merge request is not properly formatted    ###
###    the script adds the Unformatted label to the Gitlab's       ###
###    merge request, then posts in there a comment with the       ###
###    list of unformatted files asking the author to format them  ###
###                                                                ###
###  2.If there is label Unformatted and the code was formatted    ###
###    in between, then this label is removed                      ###
###                                                                ###
###  IMPORTANT: Clang-format major versions must be identical !!!  ###
###         (different versions produce different results)         ###
###                                                                ###
###      We are currently using Clang-format version 13.0.0        ###
###                                                                ###
######################################################################

main() {

  # assign basic variables
  CLANG_FORMAT_BIN="clang-format"
  CLANG_STYLE=file
  UNFORMATTED_LABEL="Unformatted"
  MERGE_REQUEST_URL="https://git.jinr.ru/api/v4/projects/$CI_MERGE_REQUEST_PROJECT_ID/merge_requests/$CI_MERGE_REQUEST_IID"
  CHANGED_FILES=""

  # do the actual work
  get_changed_files_list
  if [[ $CHANGED_FILES ]]; then
   check_formatting
   post_label_comment
  else
   echo "This merge request contains no files to format."
  fi

}


function get_changed_files_list() {

  CHANGED_FILES=$(curl -s -X GET -H "PRIVATE-TOKEN: " "$MERGE_REQUEST_URL/changes" | sed 's/,"new_path/\n&/g' | grep "new_path" \
                 | grep "\"deleted_file\":false" | cut -d '"' -f 4 | grep -E '.*\.(h|hpp|c|cpp|cxx)$' | grep -viE '.*LinkDef.h$')
  echo "Changed files to be checked for formatting:"  $CHANGED_FILES

}


function check_formatting() {

  formatted="true"
  UNFORMATTED_FILES=""
  changed_files_array=($CHANGED_FILES)

  for changed_file in "${changed_files_array[@]}"
  do
    output=$($CLANG_FORMAT_BIN -n --style=$CLANG_STYLE $changed_file |& grep -m1 warning)
    if [[ $output != "" ]]; then
      if [[ $UNFORMATTED_FILES != "" ]]; then
        UNFORMATTED_FILES="$UNFORMATTED_FILES $changed_file"
      else
        UNFORMATTED_FILES="$changed_file"
        formatted="false"
      fi
    fi
  done

  [[ $UNFORMATTED_FILES ]] && echo -e "The unformatted files are:\n $UNFORMATTED_FILES \n" || echo -e "There are no unformatted files.\n"

}


function post_label_comment() {

  FormatInfo=$(curl -s -X GET -H "PRIVATE-TOKEN: " "$MERGE_REQUEST_URL" | grep labels.*\\[.*$UNFORMATTED_LABEL.*\\])

  if [[ $formatted == false ]] && [[ -z $FormatInfo ]]; then
    curl -s -o /dev/null -S -f -X PUT -H "PRIVATE-TOKEN: $COMMENT_TOKEN" "$MERGE_REQUEST_URL?add_labels=$UNFORMATTED_LABEL"
    echo "Unformatted label added."
    author=$(curl -s -X GET -H "PRIVATE-TOKEN: " "$MERGE_REQUEST_URL" | sed 's/author/\n&/g' | grep "author" | cut -d '"' -f 11)
    comment="Dear @$author, you submitted unformatted files.<br><br>Please format them by running from your mpdroot dir:<br><br>clang-format -i **$UNFORMATTED_FILES**<br>"
    curl -s -o /dev/null -S -f -X POST -d "body=$comment" -H "PRIVATE-TOKEN: $COMMENT_TOKEN" "$MERGE_REQUEST_URL/notes"
    echo "New comment added."
  fi

  if [[ $formatted == true ]] && [[ $FormatInfo ]]; then
    curl -s -o /dev/null -S -f -X PUT -H "PRIVATE-TOKEN: $COMMENT_TOKEN" "$MERGE_REQUEST_URL?remove_labels=$UNFORMATTED_LABEL"
    echo "Unformatted label removed."
  fi

}


main "$@"
