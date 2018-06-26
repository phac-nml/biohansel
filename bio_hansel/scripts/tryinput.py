from subprocess import Popen, PIPE, call

# feht_argument_list=["./FastTree", "-gty", "nt", "f{output_directory}/core.aln", ]
#                 call()
with open('person.txt', 'w+') as file:
    call(['echo', 'hello'], stdout=file)
  

