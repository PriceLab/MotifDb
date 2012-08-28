pdf: 
	R CMD Sweave MotifDb.Rnw --pdf

view: pdf
	- kill -9 `ps aux | grep Preview | egrep -v grep | awk '{print $$2}'`
	open MotifDb.pdf

tangle:
	R CMD Stangle MotifDb.Rnw
	R --vanilla < MotifDb.R

