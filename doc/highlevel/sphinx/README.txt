(in home) 1)-3) describes setting up a virtual environment in python. Any other procedure works fine too:
1) mkdir python-environments && cd python-environments
2) virtualenv --python=python3 sphinx (sphinx here is just the name)
3) source sphinx/bin/activate
4) Install the following packages within the python environment (need to be done only once):
4.1) sudo apt-get install python3-sphinx
4.2) pip install myst-parser
4.3) pip install sphinxcontrib-mermaid
4.4) pip install sphinx-autobuild
5) Navigate to the sphinx folder (dumux/doc/highlevel/sphinx/) containing the folder build/ and source/ and execute:
   sphinx-autobuild source/ build/
6) Copy the provided address into your browser to open the html files
7) Saving the markdown files and refreshing the webpage automatically updates the html files
8) ctrl+c to stop sphinx-autobuild
9) deactivate (to close virtualenv)
