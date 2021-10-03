oMotifDb)

check: build
	(cd ..; R CMD check --no-manual --no-build-vignettes --ignore-vignettes `ls -t MotifDb_* | head -1`)

biocCheck:
	(cd ..; R CMD BiocCheck `ls -t MotifDb_* | head -1`)

unitTests: test

test:
	 for x in inst/unitTests/test_*.R; do echo ============== $$x; R -f $$x; done

site:
	R -e "devtools::build_site()"
