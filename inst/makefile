default: build install

help:
	egrep "^#" makefile | sed "s/^#//"

clean:
	- rm ../../MotifDb_*tar.gz


# --- quickbuild: no vignette
#
quickbuild:
	(cd ../..; R CMD build --no-vignettes MotifDb)


# --- build
#
build:
	(cd ../..; R CMD build --no-vignettes MotifDb)

# --- install
#
install:
	(cd ../..; R CMD install MotifDb)

# --- check
#
check: clean build install
	(cd ../..; R CMD check --no-manual --no-vignettes --no-codoc --no-examples --no-manual  MotifDb)

# --- checkfull
#
checkfull: 
	(cd ../..; R CMD build MotifDb)
	(cd ../..; R CMD check MotifDb)


# --- vanillaTest
# run all the unit tests, in a clean context
#

vanillaTest:  build install
	- rm vanillaTest.out
	R --vanilla < vanillaTest.R > vanillaTest.out 2>&1

# --- vt
# run all the unit tests, in a clean context
#

vt: vanillaTest


# --- checkvig
# check just the vignette
#

checkvig:
	(cd ../..; R CMD check --no-manual --no-codoc --no-tests --no-examples MotifDb)


# --- tangle
# extract the R code from the vignette file
#

tangle:
	(cd ../vignettes; R CMD Stangle MotifDb.Rnw)



# --- sweave
# creates MotifDb.tex, runs all embedded examples
# run this before the pdf target
#
sweave: 
	(cd ../vignettes; R CMD Sweave MotifDb.Rnw --pdf)

# --- pdf
# make and open MotifDb.pdf, the vignette 
#

pdf:  sweave
	(cd ../vignettes; open MotifDb.pdf)



