CPP	=      g++
CC      =      gcc
CFLAGS	=       -O -ffast-math     
LIBS	=         
LFLAGS	=      
VPATH=
INCLUDES = -I../../include
DEPENDS= FFT.o ran1.o rk4.o root1d.o rootNd.o spline.o lineq.o romberg.o 

default::   mathtools	

clean:
	rm -f *.o *.a  	

mathtools:   $(DEPENDS)
	ar rs ../../lib/libmathtools.a $(DEPENDS) 

FFT.o: FFT.cpp  
	$(CPP) $(INCLUDES) -c FFT.cpp $(CFLAGS) 

ran1.o: ran1.cpp  
	$(CPP) $(INCLUDES) -c ran1.cpp $(CFLAGS)

rk4.o: rk4.cpp  
	$(CPP) $(INCLUDES) -c rk4.cpp $(CFLAGS)

root1d.o: root1d.cpp  
	$(CC) $(INCLUDES) -c root1d.cpp $(CFLAGS)

rootNd.o: rootNd.c  
	$(CPP) $(INCLUDES) -c rootNd.c $(CFLAGS)

spline.o: spline.cpp  
	$(CPP) $(INCLUDES) -c spline.cpp $(CFLAGS)

lineq.o: lineq.cpp  
	$(CPP) $(INCLUDES) -c lineq.cpp $(CFLAGS)
	
romberg.o: romberg.c  
	$(CC) $(INCLUDES) -c romberg.c $(CFLAGS)
