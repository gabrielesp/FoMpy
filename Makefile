#release:
#	@echo "Uploading files to git repo"
#	git add ./*
#	git commit -m "description of version"
#	git tag -a v$(python setup.py --version) -m 'description of version'
#	git push origin --tags

pdf:	
	@echo "Generating pdf" && cd ./docs/ && $(MAKE) latexpdf

dist:
	@echo "Creating distribution files"
	python3 setup.py sdist bdist_wheel

upload:
	@echo "Uploading distribution files"
	python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*

clean:
	@echo "Removing distribution files"
	rm -rf ./build/
	rm -rf ./dist/
	rm -rf ./fompy.egg-info/
