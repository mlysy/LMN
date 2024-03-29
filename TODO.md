- [x] Uniform formatting of `@title`.  Should be: first word only capitalized, ends with ".".

- [x] `@return` lists.  These should all be of the form

    ```r
    #' @return something something a list with elements: (don't forget the ":")
    #' \describe{
    #'   \item{\code{a}}{Description of element a.}
    #'   \item{\code{b}}{Description of element b.}
    #' }
    ```

- [x] Latex code.  Typically this doesn't display properly in vim/emacs & rstudio.  So the trick is to define such commands twice, , e.g.:

    ```r
    #' \deqn{
    #' \boldsymbol{\Lambda}(\boldsymbol{Y}_t)^{\alpha - 1}
    #' }{
    #' \Lambda(Y_t)^(\alpha-1)
    #' }
    ```

    To see the different outputs, compile the PDF doc (this can be done from the terminal using `R CMD check LMN_XYZ.tar.gz`), and also look at output in Rstudio (which is a bit nicer than vim/emacs, but still very basic and it's what most people will be using)

- [x] `@description`.  These can sometimes be made more informative, if the `@title` is not enough.

- [x] Short examples for all exported functions (i.e., having the `@export` keyword).  These go in a block called `@examples`.  Also, the main doc for the package (in the file `LMN.R`) should have an example.

- [x] Package vignette.  A basic walkthrough, typically using a specific example.  I'll get this started.

- [x] Use documentation templates, i.e., `man-roxygen` and `examples`.

- [x] Pass `lmn_suff` objects to other functions, instead of entire interface.

- [x] Replace "." by "_".  

- [x] Convert roxygen documentation to Markdown formatting.

- [x] Also, consider renaming `Beta.hat` to `Bhat`.

- [x] Replace `eqn` and `deqn` roxygen tags with `preformatted`.  Actually, no.  Just write these properly.  Not so hard.

- [x] Remove `debug` arguments.

- [x] ~~Package name -- should it be `lmn`???~~

	No, I don't think this is worth doing at this point.

- [x] Use `solveV` instead of `solve` for variance matrices. 

- [ ] Translate to C++ :)

	Probably not worth the effort any time soon...

- [x] Remove `old` files (under `git` version control anyways). 

- [ ] Remove `src/DurbinLevison*` methods, and call them from **SuperGauss** instead.

	This means exposing the DL IP method in SuperGauss to call it from here...
