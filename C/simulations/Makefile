# define the C compiler to use
CC = gcc
INC = -I../headers
OPC = -g 
CFLAGS = -c -Wall -lm -O3
LDFLAGS = -lm -O3 

all: nmse_vs_snr_gaussian nmse_vs_snr_econtaminated rsnr_vs_implevel simulation_integer DCT_integer

#---------------------------------------------------#

nmse_vs_snr_gaussian: nmse_vs_snr_gaussian.o functions.o
	$(CC) $(OPC) -o $@ $^ $(LDFLAGS)

nmse_vs_snr_econtaminated: nmse_vs_snr_econtaminated.o functions.o
	$(CC) $(OPC) -o $@ $^ $(LDFLAGS)

rsnr_vs_implevel: rsnr_vs_implevel.o functions.o
	$(CC) $(OPC) -o $@ $^ $(LDFLAGS)

simulation_integer: simulation_integer.o functions.o
	$(CC) $(OPC) -o $@ $^ $(LDFLAGS)

DCT_integer: DCT_integer.o functions.o
	$(CC) $(OPC) -o $@ $^ $(LDFLAGS)

#---------------------------------------------------#

nmse_vs_snr_gaussian.o: nmse_vs_snr_gaussian.c
	$(CC) $(INC) $(CFLAGS) -c $^

nmse_vs_snr_econtaminated.o: nmse_vs_snr_econtaminated.c
	$(CC) $(INC) $(CFLAGS) -c $^

rsnr_vs_implevel.o: rsnr_vs_implevel.c
	$(CC) $(INC) $(CFLAGS) -c $^

simulation_integer.o: simulation_integer.c
	$(CC) $(INC) $(CFLAGS) -c $^

DCT_integer.o: DCT_integer.c
	$(CC) $(INC) $(CFLAGS) -c $^


#library (Don't touch)
#---------------------------------------------------#

functions.o: ../sources/functions.c 
	$(CC) $(INC) $(CFLAGS) -c $^

#---------------------------------------------------#

clean:
	rm -f *.o nmse_vs_snr_gaussian nmse_vs_snr_econtaminated rsnr_vs_implevel simulation_integer DCT_integer
