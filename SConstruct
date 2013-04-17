#-*-Python-*-

import os
import sys

compiler = ARGUMENTS.get('compiler','gcc')

if (compiler == 'gcc'):
	env = Environment(CC='gcc', CXX='g++')
elif (compiler == 'icc'): 
	env = Environment(CC='icc', CXX='icpc')
elif (compiler == 'clang'): 
	env = Environment(CC='clang', CXX='clang++')
	  
env.Append(ENV = {'PATH' : os.environ['PATH']})

env.Append(CCFLAGS='-Wall')

if (sys.platform != "darwin"):
   env.Append(LIBS='rt')


debug = ARGUMENTS.get('debug', 0)
if int(debug):
    env.Append(CCFLAGS = '-g -O0')
else:
	env.Append(CCFLAGS='-O3')
	


env.Program('test_queries', ['matrix.cpp', 'oracle_monge_matrix.cpp', 'max_value.cpp', 
							'range.cpp', 'range_query.cpp', 'envelope.cpp', 'envelope_tree.cpp', 
							'tests.cpp', 'main.cpp'])