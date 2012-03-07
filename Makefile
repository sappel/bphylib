CPP	=      g++
CC      =      gcc
CFLAGS	=       -O3 -ffast-math 
FFLAGS	=       -O  
LIBS	=         
LFLAGS	= 
VPATH   =  
INCLUDES =  -I./include 
DEPENDS = ./obj/Grid.o ./obj/poisson.o ./obj/IBS.o ./obj/target.o ./obj/landau.o ./obj/vavilov.o ./obj/Impedance.o  ./obj/Matching.o ./obj/pic.o

default::   bphylib	

clean:
	rm -f obj/*.o lib/libbphylib.a  	

bphylib:   $(DEPENDS)
	ar rs ./lib/libbphylib.a $(DEPENDS) 

obj/Grid.o: src/Grid.cpp  
	$(CPP) $(INCLUDES) -o obj/Grid.o -c src/Grid.cpp $(CFLAGS)
	
obj/poisson.o: src/poisson.cpp  
	$(CPP) $(INCLUDES) -o obj/poisson.o -c src/poisson.cpp $(CFLAGS)	
	
obj/IBS.o: src/IBS.cpp  
	$(CPP) $(INCLUDES) -o obj/IBS.o -c src/IBS.cpp $(CFLAGS)
	
obj/Impedance.o: src/Impedance.cpp  
	$(CPP) $(INCLUDES) -o obj/Impedance.o -c src/Impedance.cpp $(CFLAGS)	
	
obj/target.o: src/target.cpp  
	$(CPP) $(INCLUDES) -o obj/target.o -c src/target.cpp $(CFLAGS)
	
obj/landau.o: tools/landau/landau.f  
	gfortran $(INCLUDES) -o obj/landau.o -c tools/landau/landau.f $(FFLAGS)
	
obj/vavilov.o: tools/landau/vavilov.f  
	gfortran $(INCLUDES) -o obj/vavilov.o -c tools/landau/vavilov.f $(FFLAGS)
	
obj/Matching.o: src/Matching.cpp  
	$(CPP) $(INCLUDES) -o obj/Matching.o -c src/Matching.cpp $(CFLAGS)
	
obj/pic.o: src/pic.cpp  
		$(CPP) $(INCLUDES) -o obj/pic.o -c src/pic.cpp $(CFLAGS)
