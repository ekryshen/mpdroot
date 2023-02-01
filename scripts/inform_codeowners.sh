#!/bin/bash

# Codeowners for Gitlab v1.0
# Original version taken from GSI Darmstadt, first commit F.Uhlig
# Adapted and modified for MPDRoot @JINR Dubna by Slavomir Hnatic

#######################################################################
##################### Added functionality #############################
###                                                                 ###
### 1. The monitoring of code review process is managed by the bot. ###
###    The bot not only sets the CodeOwners label and posts         ###
###    initial comment, but also monitors code ownership changes,   ###
###    and places warning message box into description              ###
###								    ###
### 2. No git commands are used, required information is grabbed    ###
###    from Gitlab API, which removes git dependency and is faster  ###
###                                                                 ###
#######################################################################


main() {

  # Since bash only supports dictionaries from version 4 the dictionary is implemented as list. An example of dictionary entry is "/src/:roleg;hnatics".
  # Key (directory) is separated from value (list of directory owners) by a colon. Directory owners are separated by a semicolon.

  declare -a dictionary
  declare -a unique_users

  parse_codeowners_file

  MERGE_REQUEST_URL="https://git.jinr.ru/api/v4/projects/$CI_MERGE_REQUEST_PROJECT_ID/merge_requests/$CI_MERGE_REQUEST_IID"

  CHANGED_FILES=$(curl -s -X GET -H "PRIVATE-TOKEN: $COMMENT_TOKEN" "$MERGE_REQUEST_URL/changes" \
                     | sed 's/"old_path/\n&/g' | grep "old_path" | cut -d '"' -f 4)

  echo "Changed files in this merge request are:"  $CHANGED_FILES

  ALL_USERS=""
  generate_codeowners_list "$CHANGED_FILES"

  generate_comment

}

# Add a key value pair to the dictionary
function add_hash () {
  dictionary+=("$1":"$2")
}

# extract the directory (key) from the passed string
function get_key () {
  dir="${1%%:*}"
}

# extract the user list (value) from the passed string
function get_value () {
  users="${1##*:}"
}

# Parse the CODEOWNERS file and store the information of owned directories in a dictionary
# CODEOWNERS format is: entry line example /src/ @roleg @hnatics
function parse_codeowners_file () {
  while read line; do
    # exclude commented and empty lines
    if [[ $line =~ ^\# || $line =~ ^$ ]]; then
     dummy=$line
    else
      # split the line at blank characters and store the results in an array
      IFS=' ' read -ra my_array <<< "$line"
      # first entry is the directory
      directory="${my_array[0]}"
      # in case of two entries the second is the user
      if [[ ${#my_array[@]} -eq 2 ]]; then
        users="${my_array[1]}"
      # in case of more entries the entries from the second to the last are a list of users.
      # Add the users to a semicolon separated list and store the list in the user array
      else
        arr2=("${my_array[@]:1:${#my_array[@]}}")
        newUserList=""
        for i in "${arr2[@]}"; do
          newUserList="$newUserList;$i"
        done
        newUserList="${newUserList:1}"
        users="$newUserList"
      fi
      add_hash "$directory" "$users"
    fi

  done < CODEOWNERS
}

# Only the owners of the deepest directory match are informed
# Example: if CODEOWNERS contains lines  "/ @roleg" and "/scripts/ @hnatics", then only @hnatics is informed
function generate_codeowners_list () {
  CHANGED_FILES=$1
  all_users=""
  for file in $CHANGED_FILES; do
    file="/$file"
    inform_users=""
    longest_match=0
    for entry in "${dictionary[@]}"; do
      get_key $entry
      if [[ $file =~ ^$dir ]]; then
        get_value $entry
#        echo "The length of the match is ${#dir}"
        if [[ ${#dir} -gt $longest_match ]]; then
          inform_users=$users
        fi
      fi
    done
    echo "We have to inform user(s) $inform_users for file $file"
    all_users="$all_users;$inform_users"
  done
  ALL_USERS="${all_users:1}"
  # remove duplications from string
  IFS=';'
  read -ra users <<< "$ALL_USERS"
  for user in "${users[@]}"; do # access each element of array
    if [[ ! " ${unique_users[@]} " =~ " ${user} " ]]; then
      unique_users+=("$user")
    fi
 done
  OWNERS_LIST=${unique_users[@]}
  echo "We have to inform the following users about the code ownership: $OWNERS_LIST "
}

# generate current comment about the duty to review the code
# if CodeOwners label not set, then set it, post warning message into description and post comment in comments section about the duty to review
# if CodeOwners label set, then look for the last comment regarding the duty to review. Post the current comment if it differs from the last one
function generate_comment() {

 # generate current comment about the duty to review the code
  comment_start="Dear "
  comment_end="you have been identified as code owner of at least one file which was changed with this merge request. Please check the changes and approve them or request changes."

    for user in "${unique_users[@]}"; do
      comment+="$user,"
    done

  comment="$comment_start $comment $comment_end"
  echo $comment

  description_warning=$(cat << EOF

<pre>
DO NOT REMOVE THIS BOXED MESSAGE

Code in this merge request is subject to review by Code Owners.
The list of reviewers is posted in the comments section.
After all of the code is approved, please put a brief notice about it at the beginning of this message.

<i><b>NOTE: Merging this code without getting it approved by Code Owners will result in suspending your developer rights!!!</b></i>
</pre>
<br>
EOF
)


 # if codeowners label not set, then set it and place warning message box into description
 # if codeowners label already set, check for comment match
  MRInfo=$(curl -s -X GET -H "PRIVATE-TOKEN: " "$MERGE_REQUEST_URL" | grep labels.*\\[.*CodeOwners.*\\])

  if [[ -z $MRInfo ]]; then
    curl -s -o /dev/null -S -f -X PUT -H "PRIVATE-TOKEN: $COMMENT_TOKEN" "$MERGE_REQUEST_URL?add_labels=CodeOwners"
    post_comment="true"
    echo "Codeowners label added."
    # minor issue with no real effect - existing description is grabbed from Gitlab API, where it is stored in Unicode format
    # there can be some formatting discrepancy when reposting it in Gitlab UI, as Gitlab UI is using markdown format
    # i corrected only basic discrepancies and did not write full unicode --> markdown converter, as it is too much work/code and not really needed
    description_existing=$(curl -s -S -f -X GET -H "PRIVATE-TOKEN: " "$MERGE_REQUEST_URL" | sed 's/"description"/\n&/g' \
                           | sed 's/"created_at"/\n&/g' | grep -m 1 \"description\" |  cut -d '"' -f 4- | rev | cut -d '"' -f 6- | rev | sed 's/\\n/<br \/>/g' )
    description_new="$description_warning $description_existing"
    curl -s -o /dev/null -S -f -X PUT -d "description=$description_new" -H "PRIVATE-TOKEN: $COMMENT_TOKEN" "$MERGE_REQUEST_URL"
    echo "Warning message placed on top of the MR description."
  else
  last_comment_match=$(curl -s -X GET -H "PRIVATE-TOKEN: $COMMENT_TOKEN" "$MERGE_REQUEST_URL/notes" | sed 's/{"id/\n&/g' \
                         | grep "\"body\":\"$comment_start" | grep "@" | grep "$comment_end\",\"" | head -1 | grep $comment)
   if [[ -z $last_comment_match ]]; then
    post_comment="true"
    echo "List of codeowners for this merge request changed, reposting the comment in Gitlab with updated list."
   else
    echo "List of codeowners for this merge request did not change, no need to repost it in Gitlab."
   fi
  fi

 # if needed, add a comment which informs code owners about their duty to review
  if [[ -n $post_comment ]]; then
    curl -s -o /dev/null -S -f -X POST -d "body=$comment" -H "PRIVATE-TOKEN: $COMMENT_TOKEN" "$MERGE_REQUEST_URL/notes"
  fi
}

main "$@"
