import sys
import os
import shutil
import gzip
import subprocess
from multiprocessing import Pool, set_start_method
from Bio.SeqIO import index_db

parentdir = os.path.dirname(os.path.realpath(__file__))
PhyloAlndir = os.path.dirname(parentdir)
epath = os.environ['PATH']

# check if the programs exist in the path
def check_programs(progs):
	for program in progs:
		fullpath = shutil.which(program)
		if not fullpath:
			print("\nError: fail to find {}! Please install and add to it into the path!".format(program))
			sys.exit(1)

# get the path of assemble or assemble.py
def check_assemble():
	if os.path.exists(os.path.join(parentdir, 'assemble')):
		return os.path.join(parentdir, 'assemble')
	else:
		return os.path.join(parentdir, 'assemble.py')

# run the command
def runcmd(cmd, log, env=None, stdout=True, error=False):
	if env:
		# add the environment into the path
		os.environ['PATH'] = env + ':' + epath
		if stdout:
			print('PATH: +' + env)
		log.write('PATH: +' + env + "\n")
	log.write(' '.join(cmd) + "\n")
	if stdout:
		print(' '.join(cmd))
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	# synchronous output to the log file and the screen
	while p.poll() is None:
		cmdout = p.stdout.readline().decode("utf8")
		if cmdout:
			log.write(cmdout)
			if stdout:
				print(cmdout.rstrip("\n"))

	# restore to the original environment
	os.environ['PATH'] = epath

	# if error occurs, stop the script and exit
	if p.returncode != 0:
		if error:
			print("\nError in '" + ' '.join(cmd) + "'!")
			sys.exit(1)
		else:
			# if not 'error', return the state
			return False
	elif not error:
		return True

# run the function by single process with pbar
def run_sp(function, args_list, kwds={}, total=None, finish=0):
	# If total number is set, the number will be use, otherwise the length of the list will be used as total number
	if total is None:
		total = len(args_list)
	results = []

	for args in args_list:
		if kwds:
			# add the fixed parameters
			result = function(*args, **kwds)
		else:
			result = function(*args)
		results.append(result)
		# refresh the pbar
		finish += 1
		sys.stdout.write("[{}] {}{}/{} ({}%)\r".format(("+" * int((finish/total)*50)) + (" " * (50 - int((finish/total)*50))), (" " * (len(str(total))-len(str(finish)))), finish, total, '%.2f' %(finish/total*100)))
		sys.stdout.flush()

	return results

# run the function by multiprocess with pbar
def run_mp(function, args_list, cpus, kwds={}, sp_list=[]):
	multinum = len(args_list)
	total = multinum + len(sp_list)
	results = []

	if cpus == 1:
		# when using 1 cpu, use the function directly instead of multiprocess
		results.extend(run_sp(function, args_list, kwds=kwds, total=total))
	elif multinum > 0:
		# use multiprocess
		# setup multiprocess method using less memory
		try:
			set_start_method('spawn')
		except RuntimeError:
			# avoid setting repeatedly
			pass
		# setup pool
		p = Pool(cpus)
		finish = 0
		# setup results and split over cpus
		for args in args_list:
			if kwds:
				# add the fixed parameters
				results.append(p.apply_async(function, args=args, kwds=kwds))
			else:
				results.append(p.apply_async(function, args=args))
		# print the pbar
		while True:
			finish_task = sum(1 for result in results if result.ready())
			if finish_task == finish:
				continue
			# when new task finished, refresh the pbar
			finish = finish_task
			sys.stdout.write("[{}] {}{}/{} ({}%)\r".format(("+" * int((finish/total)*50)) + (" " * (50 - int((finish/total)*50))), (" " * (len(str(total))-len(str(finish)))), finish, total, '%.2f' %(finish/total*100)))
			sys.stdout.flush()
			if finish == multinum:
				break
		p.close()
		p.join()
		# extract the outputs of the function
		results = [result.get() for result in results]

	# finally manage the list of single process if set
	if sp_list:
		results.extend(run_sp(function, sp_list, kwds=kwds, finish=multinum, total=total))

	# end the pbar
	print("\n")
	return results

# read the FASTQ or FASTA file
def read_fastx(fastx, file_format='guess', select_list=None, return_iter=False, temp_idx=False):
	seqs = {}
	if file_format.endswith('_db'):
		file_format = file_format.split('_db')[0]
		if temp_idx:
			idx = 'temp.idx'
		else:
			idx = fastx + '.idx'
		try:
			db_dict = index_db(idx, fastx, file_format)
		except ValueError:
			# avoid unfinished/partial database due to last interruption
			os.remove(idx)
			db_dict = index_db(idx, fastx, file_format)
		# return the database iteration instead of directly reading the sequences
		if return_iter:
			return db_dict
		for seqid, seqinfo in db_dict.items():
			if select_list is None or seqid in select_list:
				seqs[seqid] = str(seqinfo.seq)
		db_dict.close()
		return seqs
	elif fastx.endswith('.gz'):
		reads = gzip.open(fastx, 'rt')
	else:
		reads = open(fastx)
	if file_format == 'guess':
		for line in reads:
			if line.startswith('@'):
				file_format = 'fastq'
			elif line.startswith('>'):
				file_format = 'fasta'
			break
		# recover the file iteration
		if fastx.endswith('.gz'):
			reads = gzip.open(fastx, 'rt')
		else:
			reads = open(fastx)
		print("File: '{}'\tDetected format: {}".format(fastx, file_format))
	if return_iter:
		# return the file iteration and file format instead of reading the sequences
		return reads, file_format
	seqid = None
	line_num = 0
	for line in reads:
		line = line.rstrip()
		if file_format == 'fastq' and line.startswith('@') and line_num % 4 == 0:
			seqid = line.replace('@', '', 1).replace(' ', '_')
			if select_list is not None:
				if seqid not in select_list:
					seqid = None
		elif file_format == 'fasta' and line.startswith('>'):
			arr = line.split(" ")
			seqid = arr[0].lstrip('>')
			if select_list is None or seqid in select_list:
				seqs[seqid] = ''
			else:
				seqid = None
		elif seqid:
			if file_format == 'fastq':
				seqs[seqid] = line
				seqid = None
			else:
				seqs[seqid] += line
		line_num += 1
	return seqs

