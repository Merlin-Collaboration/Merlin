# Merlin Doxygen Generation

These instructions are give for you to generate your own `doxygen`
files for `Merlin`

## Viewing the Documentation

If you wish to view pre-generated documentation, they are *temporarily*
being hosted at <http://users.ox.ac.uk/~scat6524/merlin-doxygen/html>
and a PDF copy at
<https://github.com/jfallon1997/merlin-cmake-b/releases> as well as
<http://users.ox.ac.uk/~scat6524/merlin-doxygen/latex/refman.pdf>

The documentation can be built by ensuing that BUILD_DOCUMENTATION is
enabled in cmake, and running:

    make doxygen

in the build directory. The resulting documentation will be in the
generated-docs directory.

The HTML documentation can then be viewed by opening a page (e.g.
`index.html`) in `generated-docs/html`

Man-pages and pdf-diagrams can be found and viewed in a similar manner,
under the folders `doxygen/generated-docs/man/man3` and
`doxygen/generated-docs/latex` respectively.

To generate the PDF manual, `cd` into the `latex` sub-directory and run
`make`. This will generate `refman.pdf`.

Note: to view the man pages, run the command `man
$PATH_TO_MAN_PAGES/the_page.3` or `man ,/the_page.3` *if you are in the
man page directory*.
