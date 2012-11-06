comp1:
	cp ../p1.txt .
	../gentest 20 250 11 0 250 22 255
	../../test 1 11 s 1 ll 0 0 2 > resslicot2.txt
	../../test 1 11 s 1 ll 1 0 2 > reselw2.txt
	R CMD BATCH --slave --vanilla myplots.R

