#!/usr/bin/env python
from __future__ import division, print_function
import os
import sys
import difflib
import subprocess

PY2 = sys.version_info[0] == 2
enc_arg = {"encoding":"UTF-8"}
if PY2:
	enc_arg = {}

script_dir = os.path.dirname(__file__)
uncrustify_args = "-q -l CPP -c %s/uncrustify.conf"%script_dir
uncrustify_exe = "uncrustify"

try:
	version = subprocess.check_output([uncrustify_exe, "--version"])
except:
	print("Can't find uncrustify:", uncrustify_exe)
	exit(127)

format_functions = []
def add_function(f):
	format_functions.append(f)
	return f

@add_function
def run_uncrustify(in_text):
	cmd = [uncrustify_exe] + uncrustify_args.split()
	aproc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, **enc_arg)
	as_out, as_err = aproc.communicate(input=in_text)
	if aproc.wait() != 0:
		print(uncrustify_exe, "failed to run")
		exit(1)
	return as_out

def get_whitespace(s):
	first_char = len(s) - len(s.lstrip())
	return s[:first_char], s[first_char:]

@add_function
def indent_comment_like_eclipse(in_text):
	"In block comments align leading '*' with star of '/*'"

	out_text = []
	in_block_comment = False
	for line in in_text.splitlines():
		#print(line)
		line = line.rstrip("\n")
		white, rest = get_whitespace(line)

		if rest.startswith("/*"):
			in_block_comment = True
			block_comment_firstline = True
			current_indent = white
		else:
			block_comment_firstline = False

		if "*/" in rest:
			block_comment_lastline = True
		else:
			block_comment_lastline = False

		if not in_block_comment or block_comment_firstline:
			out_text.append(line)
		elif not block_comment_lastline:
			if rest.startswith("*"):
				comment = rest.lstrip("*")
				if comment.startswith("/"): comment = " "+comment
				if comment and comment[0] not in [" ", "\t"]: comment = " "+comment
				out_text.append(current_indent + " *" + comment)
			else:
				#out_text.append(current_indent + " " + rest)
				out_text.append(line)
		else:
			out_text.append(current_indent + " " + rest)

		if block_comment_lastline:
			in_block_comment = False
	return "\n".join(out_text)+"\n"

if __name__ == "__main__":
	all_good = True
	mode = sys.argv[1]
	for fname in sys.argv[2:]:
		in_text = open(fname, **enc_arg).read()
		out_text = in_text
		#print(fname)
		for func in format_functions:
			#print(func)
			out_text = func(out_text)

		if mode == "check":
			if in_text != out_text:
				print(fname, "needs reformatting")
				all_good = False
		elif mode == "diff":
			if in_text != out_text:
				print(fname)
				sys.stdout.writelines(difflib.unified_diff(in_text.splitlines(True), out_text.splitlines(True)))
				print()
		elif mode == "fix":
			if in_text != out_text:
				open(fname, "w").write(out_text)

	if not all_good:
		exit(1)
