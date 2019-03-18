Updating ``jwst_magic`` on your machine
===================================

1. Activate your Python 3 (preferably Astroconda) environment. For installing Astroconda see: http://stsci-env.readthedocs.io/en/latest/installing_anaconda.html
2. In the directory where you last installed or updated the tools::

      $ git remote –v
3. If the result contains “JWST-FGS/Commissioning-tools.git”:

   a. You are using the old repo. We will need to clone the tools from the new repo::
   
      $ cd /location/where/you/want/the/repo
      $ git clone git@grit.stsci.edu:wfsc/tools.git
      
   b. And install the ``jwst_magic`` package::

       $ cd tools/fgs-commissioning
       $ pip install -e .
      
4. Otherwise, if the result contains “wfsc/tools.git”:

   a. You are using the tools in the new repo. Download the latest version of the repository::
   
      $ cd tools/fgs-commissioning
      $ git checkout master
      $ git pull origin master
      
   b. If your last update of the repository was before May 8th, 2018, you also need to re-install the package. This is because we added new package dependencies that might not be installed on your machine yet. (If you have updated since then though, this won’t hurt to do anyway)::

       $ pip install -e .


Running the tools
=================

Launch an IPython terminal. In IPython::

    import jwst_magic
    jwst_magic.run_tool_GUI()




