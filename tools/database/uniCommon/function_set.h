//============================================================================
// Name        : function_set.h
// Author      : Konstantin Gertsenberger (gertsen@jinr.ru)
// Description : set of common C++ functions
// Version     : 1.06
//============================================================================

#ifndef FUNCTION_SET_H
#define FUNCTION_SET_H

// C++ includes
#include <unistd.h>
#include <iomanip>
#include <iostream>
#include <cstring>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <ctime>
#include <math.h>
//#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

/* declarations */
enum BATCH_SYSTEM_NAME {
   SGE_BATCH_SYSTEM,
   TORQUE_BATCH_SYSTEM,
   SLURM_BATCH_SYSTEM
}; // 'get_batch_processor_count' function

/* OS FUNCTIONS */
// execute system command in shell (bash)
int system_command_linux(string aCommand, string &result);
// change tilda symbol on HOME in linux path
string replace_home_symbol_linux(string path);
// change $VMCWORKDIR symbol in linux path
string replace_vmc_path_linux(string path);
// get processor core count on this machine
int get_linux_processor_count();
// get maximum available processor count on Sun Grid Engine system
int get_batch_processor_count(BATCH_SYSTEM_NAME batch_system, string queue_name = "");

/* GLOBAL APPLICATION FUNCTIONS */
// get application name in linux
string get_app_name_linux();
// get aplication directory (path without file name) in linux
string get_app_dir_linux();

/* NUMBER FUNCTIONS */
// check bit in 'variable' at 'position'
#define CHECK_BIT(variable, position) ((variable) & (1ULL << (position)))

/* STRING FUNCTIONS */
// convert double number to string with a given precision
// is_fixed_point: true - {precision = number of decimal digits after point}; false - {precision = number of all decimal
// digits in the number}
string double_to_string(double number, int precision, bool is_fixed_point = true);
// convert integer number to string
string int_to_string(int number);
// convert integer (hexadecimal value) to string with hexadecimal presentation without "0x"
string int_to_hex_string(int number);
// convert string with hexadecimal presentation without "0x" to integer
int hex_string_to_int(string hex_string);
// convert string specified size in bytes to double value; "convert_to" - possible values "BKMGTP"
double byte_size_to_double(string byte_size_in_string, char convert_to = 'B');
// is string an integer number?
bool is_string_number(const string &s);
// extract first number or last number in string, only positive number by default (result number has string type)
string find_first_number(string const &str, bool isOnlyPositive = true);
string find_first_double_number(string const &str, bool isOnlyPositive = true);
string find_last_number(string const &str, bool isOnlyPositive = true);
string find_last_double_number(string const &str, bool isOnlyPositive = true);
// extract first number or last number in string, only positive number by default (result number has string type)
// beg_pos - position of the first character in the string to be assigned before the search for find_first...; returns
// position of number (string::npos, if not found) end_pos - position of the last character in the string to be assigned
// before the search for find_last...; returns position of number (string::npos, if not found)
string find_first_number(string const &str, size_t &beg_pos, size_t &end_pos, bool isOnlyPositive = true);
string find_first_double_number(string const &str, size_t &beg_pos, size_t &end_pos, bool isOnlyPositive = true);
string find_last_number(string const &str, size_t &beg_pos, size_t &end_pos, bool isOnlyPositive = true);
string find_last_double_number(string const &str, size_t &beg_pos, size_t &end_pos, bool isOnlyPositive = true);
// convert array of chars to the new lowercase array
char *convert_pchar_to_lowercase_new(char *input_char_array);
// replace string 'old_substring' by string 'new_substring' in 'text'
void replace_string_in_text(string &text, string old_substring, string new_substring);
// replace string 'old_substring' by integer 'new_subinteger' in 'text'
void replace_string_in_text(string &text, string old_substring, int new_subinteger);
// replace char 'find' in array of characters (char*) by another char 'replace'; return number of replacement
int replace_char(char *&str, char find, char replace);
// return string without leading and trailing spaces and tabs
string trim(const string &str, const string &whitespace = " \t\r");
// return string changing whitespaces and tabs by single whitespace
string reduce(const string &str, const string &fill = " ", const string &whitespace = " \t\r");
// is string ('full_str') ending with the given substring ('ending')
bool endswith(string const &full_str, string const &ending);

/*   DIR & FILE FUNCTIONS   */
// check directory exists: 0 - not exists, 1 - exists, -2 - cannot access
int check_directory_exist(const char *path);
// check and create directory if not exists: 0 - not existed before, 1 - existed, -1 - errno error, -2 - cannot access
int create_directory(const char *path);
// get file name without extension from a path
string get_file_name(string path);
// get file name with extension from path
string get_file_name_with_ext(string path);
// get directory path without last slash from file path
string get_directory_path(string file_path);

/*  TIME FUNCTIONS  */
// get current date as string
string get_current_date();
// convert string in a given format to datetime struct tm
tm convert_string_to_datetime(string str_datetime, const char *format = "%d.%m.%Y %H:%M:%S");
// convert time in double (fractional number of seconds) to a timespec
struct timespec convert_double_to_timespec(double sec);
#endif // #ifndef FUNCTION_SET_H
