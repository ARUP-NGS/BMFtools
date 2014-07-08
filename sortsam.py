def sam_sort(insam,outsam):
	#Skip header and funnel all reads into a temp file
	import random
	output=open(outsam,'w+')
	tmpname='temp{}.txt'.format(random.randint(0,200))
	tmp=open(tmpname,'w',0)
	import subprocess
	command_str=str('grep -v "@SQ\|@PG\|VN:\|@HD" {}'.format(insam))
	print(command_str)
	subprocess.call(command_str,stdout=tmp,shell=True)
	tmp.close()
	#Save the header to the outsam
	command_header = 'grep "@SQ\|@PG\|@HD" {}'.format(insam)
	subprocess.call(command_header, stdout=output, shell=True)
	#sort the reads by query name
	tmp=open(tmpname,'r')
	command_str1=str('sort -k1,1 -t " " {}'.format(tmpname))
	print(command_str1)
	subprocess.call(command_str1,stdout=output,shell=True)
	output.close()
	tmp.close()
	subprocess.call('rm {}'.format(tmpname),shell=True)
	both_cmds=command_str+"\n"+command_str1
	return(both_cmds)
