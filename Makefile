OBJ=  generate-dis.o calc_f1.o CJpdf.o F1F209.o grid_oper.o radgen.o mycernlib.o
FOR   =  -lstdc++  -lgfortran
CERNLIBS =  -L/apps/cernlib/x86_64_rhel6_4.7.2/2005/lib -lmathlib  -lpacklib 
 
semifall : $(OBJ)
	 gfortran  -o	inclusive-dis-rad  $(OBJ) $(FOR)
$(OBJ) : %.o: %.F
#	gfortran  -Dtestcj=def   -c $< -o $@  
	gfortran    -c $< -o $@

clean:
	rm -f inclusive-dis-rad   $(OBJ)

#	 gfortran  -o	inclusive-dis-rad  $(OBJ) $(CERNLIBS) $(FOR)


















