.PHONY: all style lint cov tinytex archlinux cloud build clean xz data

all: clean build

doc: README.Rmd R/*
	Rscript .dependencies_r_dev.R
	Rscript .dependencies_r_pkg.R
	Rscript -e "rmarkdown::render('README.Rmd')"
	Rscript -e "devtools::document()"

check: README.Rmd R/*
	Rscript .dependencies_r_dev.R
	Rscript .dependencies_r_pkg.R
	Rscript -e "rmarkdown::render('README.Rmd')"
	Rscript -e "devtools::document()"
	Rscript -e "devtools::check(cran = FALSE)"

cran: README.Rmd R/*
	Rscript .dependencies_r_dev.R
	Rscript .dependencies_r_pkg.R
	Rscript -e "rmarkdown::render('README.Rmd')"
	Rscript -e "devtools::document()"
	Rscript -e "devtools::check()"

style:
	Rscript -e "styler::style_dir(exclude_dirs = '.library', filetype = c('.R', '.Rmd'))"

lint:
	Rscript -e "lintr::lint_dir('R')"
	Rscript -e "lintr::lint_dir('tests')"

cov:
	Rscript -e "covr::package_coverage()"

tinytex:
	Rscript .dependencies_tinytex.R

archlinux: README.Rmd R/* vignettes/*
	bash .dependencies_archlinux_dev.sh
	bash .dependencies_archlinux_pkg.sh
	Rscript .dependencies_r_dev.R
	Rscript .dependencies_r_pkg.R

cloud: README.Rmd R/* vignettes/*
	Rscript .dependencies_r_dev.R
	Rscript .dependencies_r_pkg.R

build:
	Rscript .dependencies_r_dev.R
	Rscript .dependencies_r_pkg.R
	Rscript -e "styler::style_dir(exclude_dirs = '.library', filetype = c('.R', '.Rmd'))"
	Rscript -e "rmarkdown::render('README.Rmd')"
	Rscript -e "devtools::document()"
	Rscript -e "devtools::check(cran = FALSE)"
	Rscript -e "pkgdown::build_site()"
	Rscript -e "devtools::build(path = '.')"

xz:
	Rscript -e "tools::resaveRdaFiles(paths = 'data', compress = 'xz')"

data:
	Rscript -e "tools::resaveRdaFiles(paths = 'data')"

clean:
	@rm -rf README.html
	@rm -rf README.md
	@rm -rf docs/*
	@rm -rf man/*
	@rm -rf NAMESPACE
