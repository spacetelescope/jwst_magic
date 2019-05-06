Python Code Style Guide for `jwst_magic`
==============================================

This document serves as a style guide for all `jwst_magic` software development.  Any requested contribution to the `jwst_magic` code repository should be checked against this guide, and any violation of the guide should be fixed before the code is committed to
the `master` branch.

Prerequisite Reading
--------------------

It is assumed that the reader of this style guide has read and is familiar with the following:

- The [PEP8 Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/)
- The [PEP257 Docstring Conventions Style Guide](https://www.python.org/dev/peps/pep-0257/)
- The [`numpydoc` docstring convention](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt)


Workflow
--------

All software development for the `jwst_magic` project should follow a continuous integration workflow, described in the [`git` & GitLab Workflow for Contributing](./style_guide/git_workflow.md).  Before committing any code changes, use a linter such as `pylint` and with your editor to check against PEP8, PEP257, and `numpydoc` docstring standards. Also check that your code is conforming to this style guide.


Security
--------

The following items should never be committed in the `jwst_magic` source code or GitHub issues/pull requests:

- Account credentials of any kind (e.g. database usernames and passwords)
- Internal directory structures or filepaths
- Machine names
- Proprietary data

Additionally, developers of this project should be mindful of application security risks, and should adhere to the [OWASP Top 10](https://www.owasp.org/images/7/72/OWASP_Top_10-2017_%28en%29.pdf.pdf) as best possible.


`jwst_magic`-Specific Code Standards
------------------------------

`jwst_magic` code shall adhere to the `PEP8` conventions save for the following exceptions:

 - Lines of code need not to be restricted to 79 characters.  However, it is encouraged to break up obnoxiously long lines into several lines if it benefits the overall readability of the code


`jwst_magic`-Specific Documentation Standards
---------------------------------------

`jwst_magic` code shall adhere to the `PEP257` and `numpydoc` conventions.  The following are further recommendations:

- Each function/method should have at minimum a description, `Parameters` (if necessary), and `Returns` (if necessary) sections


Attribution
------------
This workflow is adapted from the [`spacetelescope` `jwql` Python Code Style Guide for `jwql`](https://github.com/spacetelescope/jwql/blob/master/style_guide/style_guide.md).
