FC = ifort 
FFLAGS +=  -w -O3

 
AR = ar r  
LINKLIB = ld -shared  
LIBDIR = ./
 
NSOBJECTS = utils.o utils1.o priors.o kmeans_clstr.o xmeans_clstr.o posterior.o nested.o

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^ 

%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $^ 

 
all: libnest3.a  
 
libnest3.so: $(NSOBJECTS) 
	$(LINKLIB) -o $(LIBS) $@ $^ 
 
libnest3.a: $(NSOBJECTS) 
	$(AR) $@ $^ 
 
clean: 
	-rm $(LIBDIR)/libnest3.*  *.o *.mod
	
