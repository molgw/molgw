#!/usr/bin/env python3

import os, sys, glob, shutil
import filecmp
from pathlib import Path

MAX_LINE_LENGTH = 132
INDENT_UNIT = 2

if len(sys.argv) < 2:
    print('Please provide a file or a folder')
    sys.exit(1)
else:
    path = sys.argv[-1]

if ".f90" in path:
    ffiles = [path]
else:
    ffiles = glob.glob(path+'/*f90')

ffiles = [str(Path(ffile).resolve()) for ffile in ffiles ]

print(ffiles)
new_directory = ffiles[0].rsplit("src/")[0] + "src_formatted/"
print(f"Target directory: {new_directory}")
if os.path.exists(new_directory):
    shutil.rmtree(new_directory)
os.makedirs(new_directory)
#Path(new_directory).mkdir(parents=True, exist_ok=True)

increase_indent = ["program", "module", "subroutine", "function", "do", "if", "select", "block", "interface", "type", "else", "case", "contains"]
decrease_indent = ["end", "enddo", "endif", "endblock", "endinterface", "endselect", "endtype", "else", "case", "contains"]

def line_length(string):
    return len(line.split("!")[0])

def count_indent(string):
    return len(string) - len(string.lstrip())

def get_first_word(string):
    tmp_string = string.lstrip()
    if len(tmp_string) == 0:
        word = ""
    elif tmp_string[0] in ["#", "!", "&"]:
        word = tmp_string[0]
    else:
        word = tmp_string.split(" ")[0].split("(")[0].rstrip()
    return word


for ffile in ffiles:
    short_file_name = ffile.split("/")[-1]
    print(f"Analyze file: {short_file_name}")

    with open(ffile,"r") as fo:
        old_file_string = [line.rstrip() + line[-1]  for line in fo]

    new_file_string = []

    target_indent = 0
    continuation_line = False

    for i, line in enumerate(old_file_string):
        iline = i + 1
        #
        # Check line length
        if line_length(line) > MAX_LINE_LENGTH:
            print(f"Line {iline:05d}: too long")
            print(line.rstrip())
        #
        # Check for tabs (not allowed!)
        if "\t" in line:
            print(f"Line {iline:05d}: tab found")
            print(line.rstrip())

        first_word = get_first_word(line)
        #print(iline,"__"+first_word+"__")

        if first_word == "&":
            print(f"Line {iline:05d}: continuation sign at the beginning deprecated")
            print(line.rstrip())

        #
        # Check indentation
        #

        # Skip empty line or continuated line
        if len(first_word) == 0 or continuation_line:
            #print(f"{iline} is skipped")
            continuation_line = False
            new_file_string.append(line)

        elif first_word == "#":
            # Treat preprocessor line
            if count_indent(line) != 0:
                print(f"Line {iline:05d}:  {count_indent(line)}   0"  )
                print(line.rstrip())
            new_file_string.append(line)
        else:
            # Regular line

            if first_word in decrease_indent:
                target_indent -= INDENT_UNIT
                target_indent = max(target_indent,0)
                
            if count_indent(line) != target_indent:
                print(f"Line {iline:05d}: fix indentation found {count_indent(line)} instead of {target_indent}"  )
                print(line.rstrip())
                new_file_string.append(" "*target_indent + line.lstrip())
            else:
                new_file_string.append(line)

            if first_word in increase_indent:
                target_indent += INDENT_UNIT
                #
                # Then cancel the increase for some ambiguous cases
                #
                # There are two meanings for contains
                if first_word == "contains" and target_indent == INDENT_UNIT:
                    target_indent -= INDENT_UNIT
                # There are two meanings for type
                if first_word == "type" and "::" in line:
                    target_indent -= INDENT_UNIT

            # Treat the special case of "if" without "then"
            if first_word == "if" and not "then" in line:
                j = i
                while "&" in old_file_string[j]:
                    j += 1
                if not "then" in old_file_string[j]:
                    target_indent -= INDENT_UNIT

        if "&" in line:
            continuation_line = True

    new_file = new_directory + short_file_name
    with open(new_file,"w") as ff:
        for line in new_file_string:
            ff.write(line)

    identical = filecmp.cmp(new_file,ffile)
    if identical:
        print(f"   {ffile} is left unchanged")
        os.remove(new_file)




