
NCFCODES := $(wildcard *.ncf)
SIMBINARIES := $(patsubst %.ncf,%_simulate,$(NCFCODES))
SIMCODES := $(patsubst %.ncf, %_simulate.c,$(NCFCODES))

all: $(SIMBINARIES)

.SECONDARY: $(SIMCODES)

%_simulate.c: %.ncf neurowriter.py
	python neurowriter.py $< 
%_simulate: %_simulate.c
	gcc -pedantic -Wall -O3 -o $@ $<


#target: dependencies
#<tab>	commands to compile

#fhn_simulate: fhn_simulate.c 
#	
#izhi_simulate: izhi_simulate.c
#	gcc -o izhi_simulate izhi_simulate.c -Wall

echo:
	@echo $(NCFCODES)
	@echo $(SIMCODES)
	@echo $(SIMBINARIES)
clean:
	rm -f $(SIMCODES) $(SIMBINARIES)
