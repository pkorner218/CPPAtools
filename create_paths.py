import os
import logging
import sys

#- mains  -lib
#        - s4pred
#        - validated_fasta/
#        - chem/
#        - dom/
#        - FS_output/   -metadata


cwd = os.getcwd()


def mkdir_exc(filepathname):

	try:
		os.mkdir(filepathname)
		#logging.info(f" '{filepathname}' created successfully.")
	except FileExistsError:
		logging.info(f" '{filepathname}' already exists.")
		logging.warning("Data cannot be overwritten ! Please ensure you start in a new directory or use only the separate steps parts you require.")
		sys.exit(1)
	except PermissionError:
		logging.info(f"Permission denied: Unable to create '{filepathname}'.")
		logging.warning("Pipeline can not run without permissions. Please assure permission policy is correct.")
		sys.exit(1)
	except Exception as e:
		logging.info(f"An unknown error occurred: {e}")
		sys.exit(1)

def main():

	logging.basicConfig(level = logging.INFO , format='{asctime} - {levelname} - {message}',filemode='w', style = "{", datefmt = "[%d-%m-%Y] [%H:%M] ")

	validatedfilepath = cwd + "/validated_fasta/"
	mkdir_exc(validatedfilepath)

	aberrfilepath = cwd + "/aberr_fasta/"
	mkdir_exc(aberrfilepath)

	chemfilepath = cwd + "/chem/"
	mkdir_exc(chemfilepath)

	domfilepath = cwd + "/dom/"
	mkdir_exc(domfilepath)

	secstrucpath = cwd + "/sec_struc/"
	mkdir_exc(secstrucfilepath)

main()
