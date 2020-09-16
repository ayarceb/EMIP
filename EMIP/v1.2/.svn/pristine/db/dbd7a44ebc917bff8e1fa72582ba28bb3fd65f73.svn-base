#
# Makefile for main directory.
#

help:
	@echo " "
	@echo " Usage:"
	@echo "   make docu           # build documentation"
	@echo "   make clean          # remove temporary files"
	@echo "   make clean-all      # remove temporary files and generated documentation"
	@echo " "

docu:
	@echo " "
	(cd doc; make html)
	@echo " "
	@echo "Browse through documentation:"
	@echo "   file://$(PWD)/doc/build/html/index.html"
	@echo " "

clean: clean-pyc clean-osx clean-docu
	rm -f *.jb *.out *.err

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} \;
	rm -f -r ./py/__pycache__

clean-osx:
	find . -name '.DS_Store' -exec rm -f {} \;

clean-docu:
	rm -f -r doc/build/doctrees

clean-all: clean
	rm -f -r doc/build/html

