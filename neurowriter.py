#!/usr/bin/python
""" neurowriter """
"""
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

import sys,os
class NeuroIntegrator():
	def __init__(self, name,Symbols,Equations, Params,Settings,OutputFile):
		self.name = name
		self.sym = Symbols
		self.eqs = Equations
		self.params = Params
		self.settings = Settings
		self.output = OutputFile
	def print_1step(self):
		state = ", ".join(["double %s_start, double *%s_new" % (v, v) for (k,v) in self.sym.items()])
		sysparams = ", ".join(["double %s"%s for s in self.params.split(", ")])
		stdinparams= ", ".join(["double %s"%s for s in self.settings['stdinparams'].split(',')])

		self.output.write("void %s_1step(%s, \n\t\tdouble shared, unsigned long *i_seed, %s, %s) {\n"%(self.name,state,sysparams,stdinparams))
		self.output.write("\t/* init state */\n")
		self.output.write("\tdouble %s;\n" % ", ".join([s.upper() for s in self.sym.values()]))

                for var in self.settings['staticparams'].split(','):
			eq = self.settings['static_%s'%var.strip()]
			for (K,stateVar) in self.sym.items():
				VAR = self.sym[K].strip()
                                VAR_start = "%s_start"%VAR
                                eq = eq.replace("@%s"%VAR,VAR_start)
                        self.output.write("\tdouble %s = %s;\n" % (var.strip(),eq))
		#for (k,v) in self.sym.items():
		#	print "\tdouble %s;"%(v.upper())
		self.output.write("\t/* init rand */\n")
		self.output.write("\tdouble unique_stoch = gasdev();\n")
        # self.output.write("\t//printf(\"%f %f\\n\", shared, unique_stoch);\n")
		#print "\tlong i_seed = (unsigned int)time(NULL);"
		#print "\tdouble RANDN = 0.00;"
		for (k,v) in self.eqs.items():
			var = self.sym[k].strip()
			var_start = "%s_start"%var
			var_new = "%s_new"%var
			var_tmp = var.upper()
			for (K,V) in self.sym.items():
				VAR = self.sym[K].strip()
				VAR_start = "%s_start"%VAR
				v = v.replace("@%s"%VAR,VAR_start)
			str =  "\t%s = %s + (dt*(%s)"%(var_tmp, var_start,v)
			self.output.write("/* output: k = [%s] */\n"%k)
			if int(k)==1:
				str+=" + sqrt(dt)*sigma*(sqrt(sharedInput)*shared + sqrt(1-sharedInput)*unique_stoch)"
			if self.settings['scale_eq%s'%(k)] is not "":
				scale = "/%s"%self.settings['scale_eq%s'%(k)]
			else:
				scale = ""
			self.output.write("%s)%s;\n"%(str,scale))
			self.output.write("\t*%s = %s;\n" % (var_new, var_tmp))
		self.output.write("}\n")
	#def print_aux_functions(self):
	#	aux_Functions = [(k,v) for (k,v) in self.settings.items() if k.startswith('auxillary_function')]
	#	for (entry, function) in aux_Functions:
	#		function_name = entry[19:]
	#		print "double %s(
	def print_euler_loop(self):
		self.output.write( """\n
int main(int argc, char** argv){\n""")
        #double a=-0.7; double b=0.8;"""
		self.output.write("\n".join(["\tdouble %s_start[2], %s_new[2];"%(s,s) for (i,s) in self.sym.items()]))
		self.output.write("""\n\n
        /* double dt=0.01,tMax = 500; */\n
        /* double mu=.47; double sigma=0.10;double sharedInput = 0.50; */\n
        int i,j;\n
	double shared;\n""")
	
		self.output.write("\tif (argc < %d) { \n"%(len(self.params.split(','))+1))
		self.output.write("\t\tprintf(\"oops, not enough args\\n\");\n")
		self.output.write("\t\tprintf(\"usage: %%s %s i_seed\\n\",argv[0]);\n"% (" ".join(self.params.split(','))))
		self.output.write("\t\tprintf(\"stdin vars: %s\\n\");\n" % (" ".join(self.settings['stdinparams'].split(','))))
		self.output.write("\t\tprintf(\"configured actions: %s\\n\");\n" % (" ".join(self.settings['actions'].split(','))))
		self.output.write("\t\treturn 1;\n")
		self.output.write("\t}\n")
		self.output.write("\tchar *endp;\n")
		cmdParamCounter = 1
		for var in self.params.split(','):
			self.output.write("\tdouble %s = strtod(argv[%s],&endp);\n"%(var,cmdParamCounter))
			cmdParamCounter += 1
		#print "\tunsigned long i_seed = atoi(argv[%i]);" % (cmdParamCounter)
		self.output.write("\tunsigned long i_seed= ((unsigned long)time(NULL));\n")
		self.output.write("\tif (atoi(argv[%i]) > 0)\n"%(cmdParamCounter))
		self.output.write("\t\ti_seed = atoi(argv[%i]);\n"%(cmdParamCounter))

		self.output.write("\tinit_genrand(i_seed);\n")
		# read in stdin variables for things
		stdinVarCount = 0
		argStr = []
		for var in self.settings['stdinparams'].split(','):
			self.output.write("\tdouble %s;\n"%(var.strip()))
			argStr = argStr + ["&%s"%var.strip()]
			stdinVarCount += 1
		readStr = " ".join(["%lf"]*stdinVarCount)
		self.output.write( "\tscanf(\"%s\",%s);\n"%(readStr,", ".join(argStr)))
#		for var in self.settings['staticparams'].split(','):
#			print "\tdouble %s = %s;" % (var.strip(), self.settings['static_%s'%var.strip()])

		self.output.write("/* initial conditions */\n")
		self.output.write("\t for (j=0; j<2; j++) {\n")
		inits = lambda stage: ["\t\t%s_%s[j]=%s;"%(s, stage, self.settings['ic_var%d'%(int(i))]) for (i,s) in self.sym.items()]
		self.output.write("\n".join(inits('start')+ inits('new')))
		self.output.write("\n\t}\n")

		self.output.write("\tfor (i=0; i*dt <= tMax; i++) {\n")
                self.output.write("\t\tfor (j=0; j<2; j++) {\n")
		self.output.write("\t\t\t /* Optional Spike Detection*/\n")
		self.output.write("/* settings = %s */\n" % self.settings)

		if 'printall' in self.settings['actions']:
			#for (k,v) in self.sym.items():
			vars = lambda T: ["%s_%s[j]"%(statevar,T) for i,statevar in self.sym.items()]
			stateList = vars('new') + vars('start')
			format = "%3.15f "*len(stateList)
			state = ",".join(stateList)
			self.output.write("\t\t\tprintf(\"%%d %s\\n\",j,%s);\n"%(format,state))
		if 'printspiketimes' in self.settings['actions']:
			if self.settings['spike_trigger'] == "threshold":
				thresh = self.settings['spike_threshold']
				var = self.settings['spike_var'][1:]
				self.output.write("\t\t\tif ((%s_new[j]>%s)&&(%s_start[j]<=%s)) {\n" % (var,thresh,var,thresh))
				self.output.write("\t\t\t\tprintf(\"%d %f\\n\",j+1,dt*i);\n")
				self.output.write("\t\t\t}\n")


		if 'reset' in self.settings['actions']:
			self.output.write("/* do reset logic here! */\n")
			self.output.write("/* self.sym = %s */\n" % (self.sym))
			self.output.write("\t\tif (%s_new[j] > %s) {\n" % (self.settings['reset_variable'][1:], self.settings['reset_threshold']))
			for (k,v) in self.sym.items():
						# logic for 1 v: set v -65
						# logic for 2 u: add u 100
				style = self.settings['reset_var%s_style'%k]
				val = self.settings['reset_var%s_val'%k]
				self.output.write("\t\t\t /* logic for %s %s: %s %s %s */\n" % (k,v,style, v,val))
				if style == "add":
					self.output.write("\t\t\t%s_new[j] = %s_new[j] + %s;\n"%(v,v,val))
				elif style == "set":
					self.output.write("\t\t\t%s_new[j] = %s;\n"%(v,val))
				else:
					print "weird style... bailing"
					sys.exit(15)
			self.output.write("\t\t}\n")
#			elif self.settings['spike'] = "threshold":
#				print "\t\t\t\tif (v_new[j]>%s) {"
#				print "\t\t\t\t\t
#				print "\t\t\t\t}"
		self.output.write("\n".join(["\t\t\t%s_start[j] = %s_new[j];"%(s,s) for (i,s) in self.sym.items()]))
		self.output.write("\n\t\t}\n")
		self.output.write("\tshared = gasdev();\n")
                self.output.write("\t\tfor (j=0; j<2; j++) {\n")
		
		state = ", ".join(["%s_start[j], &%s_new[j]"%(s,s) for (i,s) in self.sym.items()])
		params = ", ".join(self.params.split(", "))
		stdin = ", ".join(self.settings['stdinparams'].split(', '))

		self.output.write("\t\t\t%s_1step(%s,\n\t\t\t\tshared, &i_seed, %s, %s);\n" % (self.name,state,params,stdin))
		self.output.write("""\n
                }\n
        }\n
        return 0;\n
}\n
""")
	def print_basics(self):
		self.output.write("""#include <stdio.h>\n
#include <time.h>\n
#include <stdlib.h>\n
#include <math.h>\n
double gasdev();\n
void init_genrand(unsigned long s);\n
""")


def main():
	from sys import argv
	f = open(argv[1]) # ncf = neuro comp file
	SymbolTable = {}
	Equations = {}
	Settings = {}
	#params = []
	Name = ""
	for line in f.readlines():
		if line[0] == ";" or  line.find("=") < 0:
			continue
		(lhs, rhs) = line.strip().split("=")
		(lhs,rhs) = (lhs.strip(), rhs.strip())
		if lhs.startswith('var'):
			index = lhs[3:]
			SymbolTable[index] = rhs
		elif lhs.startswith('eq'):
			index = lhs[2:]
			Equations[index] = rhs
			Settings['scale_eq%s'%(index)] = ""
		elif lhs.startswith('cmdparams'):
			params = rhs
		elif lhs.startswith('name'):
			Name = rhs
		else:
			#print "unspecified command: %s"% lhs
			Settings[lhs.strip()] = rhs.strip()
	outputFile = open(os.path.splitext(argv[1])[0]+'_simulate.c','w')
	thing = NeuroIntegrator(Name, SymbolTable, Equations, params, Settings, outputFile)
	thing.print_basics()
	thing.print_1step()
	thing.print_euler_loop()
	f = open('boilerplate');
	for line in f.readlines():
		outputFile.write(line)
	
if __name__ == "__main__":
	main()
