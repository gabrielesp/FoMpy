git:
	@echo "Uploading files to git repo"
	rm -rf ./fompy/__pycache__/
	git add ./*
	echo -n'Write commit description'
	read newdescription
	git commit -m "$newdescription"
	git push origin master

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
